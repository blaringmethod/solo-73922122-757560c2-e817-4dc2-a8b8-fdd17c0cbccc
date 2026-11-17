"""Index (barcode) verification for paired FASTQ task results.

Streaming dual-index demultiplexing audit for a *succeeded* pairing task.
For every paired record the first eight bases of the R1 and R2 sequences are
read as the observed indices; each end tolerates at most one mismatch against
the sample sheet. A pair is assigned to the single sample whose combined
R1+R2 mismatch count is smallest; tied optima are *ambiguous*. When no sample
is legal on both ends, a pair whose ends each uniquely point at two different
samples is a suspected sample crossover, everything else is unassigned.

The result files are read through fixed-size chunks with a small line scanner,
so records with very long sequences never have to fit in memory and the whole
process stays inside the 256 MiB application cap. Files are only opened for
reading: a failed check never mutates the original results and never returns
partial statistics.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Optional

from . import db, worker
from .fastq import FastqFormatError, extract_read_id

INDEX_LEN = 8
"""Observed index prefix length in bases."""

MAX_MISMATCH_PER_END = 1
"""At most one mismatch is allowed on either end."""

MIN_SAMPLES = 2
MAX_SAMPLES = 32

# A succeeded task flips its status immediately before the result temp files
# are renamed into place, so a check arriving in that narrow window waits
# briefly instead of reporting a spurious missing-result error.
PUBLISH_SETTLE_RETRIES = 25
PUBLISH_SETTLE_INTERVAL_S = 0.05

READ_CHUNK = 64 * 1024
"""Fixed read size; only a few chunks are live in memory at any moment."""

LINE_CAP_BYTES = 1 * 1024 * 1024
"""Header/plus lines are structurally short; reject runaway lines early."""

# 2-bit encoding of the four DNA bases.
_BASE_CODE = {ord("A"): 0, ord("C"): 1, ord("G"): 2, ord("T"): 3}
_MASK_LOW_BITS = 0x5555  # low bit of every 2-bit group, for 8 index positions


class IndexCheckError(Exception):
    """Whole-request failure of an index check (no partial statistics)."""

    def __init__(self, status_code: int, detail: str) -> None:
        self.status_code = status_code
        self.detail = detail
        super().__init__(detail)


class SampleSpec:
    """One sample sheet row with pre-encoded indices."""

    __slots__ = ("name", "r1_raw", "r2_raw", "r1_code", "r2_code")

    def __init__(self, name: str, r1_index: str, r2_index: str) -> None:
        self.name = name
        self.r1_raw = r1_index.encode("ascii")
        self.r2_raw = r2_index.encode("ascii")
        self.r1_code = _encode_index(self.r1_raw)
        self.r2_code = _encode_index(self.r2_raw)


def _encode_index(raw: bytes) -> int:
    """Encode an already validated 8-char ACGT index as a 16-bit integer."""
    code = 0
    for byte in raw:
        code = (code << 2) | _BASE_CODE[byte]
    return code


def _encode_prefix(prefix: bytes) -> Optional[int]:
    """Encode an observed prefix; return ``None`` if it contains non-ACGT."""
    code = 0
    for byte in prefix:
        value = _BASE_CODE.get(byte)
        if value is None:
            return None
        code = (code << 2) | value
    return code


def _hamming_bytes(left: bytes, right: bytes) -> int:
    return sum(1 for a, b in zip(left, right) if a != b)


def _candidates(
    prefix: bytes, codes: list[int], raw_indices: list[bytes]
) -> list[tuple[int, int]]:
    """Return ``(sample_index, distance)`` pairs with distance <= 1.

    A prefix shorter than eight bases (or an empty one) matches no sample.
    Non-ACGT characters fall back to a byte-wise Hamming distance.
    """
    if len(prefix) != INDEX_LEN:
        return []
    encoded = _encode_prefix(prefix)
    if encoded is None:
        return [
            (i, d)
            for i, raw in enumerate(raw_indices)
            if (d := _hamming_bytes(prefix, raw)) <= MAX_MISMATCH_PER_END
        ]
    out: list[tuple[int, int]] = []
    for i, sample_code in enumerate(codes):
        diff = encoded ^ sample_code
        # Collapse each 2-bit group to one flag bit, then count the groups.
        distance = ((diff | (diff >> 1)) & _MASK_LOW_BITS).bit_count()
        if distance <= MAX_MISMATCH_PER_END:
            out.append((i, distance))
    return out


def classify_pair(
    prefix_r1: bytes, prefix_r2: bytes, samples: list[SampleSpec]
) -> tuple[str, Optional[int]]:
    """Classify one read pair.

    Returns ``("sample", index)`` for a unique best assignment,
    ``("ambiguous", None)`` for a tied optimum, ``("cross_sample", None)`` for
    a suspected crossover and ``("unassigned", None)`` otherwise.
    """
    c1 = _candidates(
        prefix_r1, [s.r1_code for s in samples], [s.r1_raw for s in samples]
    )
    c2 = _candidates(
        prefix_r2, [s.r2_code for s in samples], [s.r2_raw for s in samples]
    )

    r2_dist = dict(c2)
    legal = [
        (d1 + r2_dist[i], i)
        for i, d1 in c1
        if i in r2_dist  # both ends must point at the SAME sample
    ]
    if legal:
        best_total = min(total for total, _ in legal)
        winners = [i for total, i in legal if total == best_total]
        if len(winners) > 1:
            return "ambiguous", None
        return "sample", winners[0]

    # No legal combined assignment: crossover only when each end has exactly
    # one candidate and the two candidates are different samples.
    if len(c1) == 1 and len(c2) == 1 and c1[0][0] != c2[0][0]:
        return "cross_sample", None
    return "unassigned", None


class _ChunkReader:
    """Line scanner over a binary file with a hard bound on live bytes.

    ``read_line`` never holds more than a couple of READ_CHUNK buffers: long
    sequence/quality lines are measured and only their first ``keep_prefix``
    content bytes are retained. A CR immediately before an LF (the CRLF case,
    even when the split lands between chunks) is excluded from the content;
    any other CR, including a lone CR at EOF, is ordinary content.
    """

    def __init__(self, fh, chunk_size: int = READ_CHUNK) -> None:
        self._fh = fh
        self._chunk_size = chunk_size
        self._buf = b""

    def read_line(
        self, *, cap: Optional[int] = None, keep_prefix: int = 0
    ) -> tuple[bytes, int, bool]:
        """Read through the next LF (or EOF).

        Returns ``(kept_bytes, content_length, ended_at_eof)``. With ``cap``
        the full line content is retained and a line past the cap raises
        ``FastqFormatError``; with ``keep_prefix`` only that many leading
        content bytes are retained.
        """
        parts: Optional[list[bytes]] = [] if cap is not None else None
        prefix = bytearray()
        total = 0

        def add(body: bytes) -> None:
            nonlocal total, prefix
            if not body:
                return
            total += len(body)
            if parts is not None:
                parts.append(body)
                if total > cap:
                    raise FastqFormatError(
                        f"line exceeds the {cap} byte structural limit"
                    )
            if keep_prefix and len(prefix) < keep_prefix:
                prefix += body[: keep_prefix - len(prefix)]

        while True:
            if not self._buf:
                self._buf = self._fh.read(self._chunk_size)
                if not self._buf:
                    if parts is not None:
                        return b"".join(parts), total, True
                    return bytes(prefix), total, True

            newline = self._buf.find(b"\n")
            if newline >= 0:
                body = self._buf[:newline]
                self._buf = self._buf[newline + 1 :]
                if body.endswith(b"\r"):
                    body = body[:-1]
                add(body)
                if parts is not None:
                    return b"".join(parts), total, False
                return bytes(prefix), total, False

            # No LF in the buffer. Hold back a trailing CR: if the next read
            # starts with LF it belongs to a CRLF; otherwise it is content.
            if self._buf.endswith(b"\r"):
                add(self._buf[:-1])
                self._buf = b"\r"
                nxt = self._fh.read(self._chunk_size)
                if not nxt:
                    # No LF follows: the parked CR is ordinary content.
                    add(b"\r")
                    if parts is not None:
                        return b"".join(parts), total, True
                    return bytes(prefix), total, True
                self._buf = b"\r" + nxt
            else:
                add(self._buf)
                self._buf = b""


def iter_result_records(path: Path, side_label: str):
    """Yield ``(read_id, index_prefix)`` per record, streaming every line.

    Applies the original four-line FASTQ contract: header starting with '@',
    third line starting with '+', non-empty sequence/quality and equal lengths.
    Raises ``FastqFormatError`` on any structural violation and ``OSError``
    on storage faults; the file is opened read-only and never modified.
    """
    with path.open("rb") as fh:
        reader = _ChunkReader(fh)
        record_no = 0
        while True:
            header, _, header_eof = reader.read_line(cap=LINE_CAP_BYTES)
            if header_eof and not header:
                if record_no == 0:
                    raise FastqFormatError(
                        f"{side_label}: result file contains no records"
                    )
                return
            record_no += 1
            if not header.startswith(b"@"):
                raise FastqFormatError(
                    f"{side_label}: record {record_no} header does not start "
                    "with '@'"
                )
            read_id = extract_read_id(header)

            prefix, seq_len, seq_eof = reader.read_line(keep_prefix=INDEX_LEN)
            if seq_eof and seq_len == 0:
                raise FastqFormatError(
                    f"{side_label}: record {record_no} is truncated "
                    "(expected four lines)"
                )
            if seq_len == 0:
                raise FastqFormatError(
                    f"{side_label}: record {record_no} has an empty sequence"
                )

            plus, _, plus_eof = reader.read_line(cap=LINE_CAP_BYTES)
            if (plus_eof and not plus) or not plus.startswith(b"+"):
                raise FastqFormatError(
                    f"{side_label}: record {record_no} third line does not "
                    "start with '+'"
                )

            _, qual_len, qual_eof = reader.read_line()
            if qual_eof and qual_len == 0:
                raise FastqFormatError(
                    f"{side_label}: record {record_no} is truncated "
                    "(missing quality line)"
                )
            if qual_len == 0:
                raise FastqFormatError(
                    f"{side_label}: record {record_no} has an empty quality "
                    "string"
                )
            if seq_len != qual_len:
                raise FastqFormatError(
                    f"{side_label}: record {record_no} sequence length "
                    f"{seq_len} != quality length {qual_len}"
                )

            yield read_id, prefix


def _resolve_result_paths(task_id: str, row) -> tuple[Path, Path]:
    base = worker.task_dir(task_id).resolve()
    rels = (row["result_r1"], row["result_r2"])
    if not all(rels):
        raise IndexCheckError(
            422, "task result is missing a published result file path"
        )
    resolved = [(base / rel).resolve() for rel in rels]
    # Guard against traversal even though the rel values are server-generated.
    for path in resolved:
        if base not in path.parents:
            raise IndexCheckError(422, "result path escapes the task directory")
    # Tolerate the small window between the 'succeeded' status flip and the
    # atomic rename that publishes the result files.
    for _ in range(PUBLISH_SETTLE_RETRIES):
        if all(path.is_file() for path in resolved):
            return resolved[0], resolved[1]
        time.sleep(PUBLISH_SETTLE_INTERVAL_S)
    raise IndexCheckError(
        422, f"result file is missing or unreadable: {rels[0] if not resolved[0].is_file() else rels[1]}"
    )


def run_index_check(task_id: str, samples: list[SampleSpec]) -> dict:
    """Stream-audit one succeeded task and return the full statistics dict.

    Any missing result, storage fault or corrupt record aborts the whole call
    with ``IndexCheckError``/``FastqFormatError``; counters accumulated so far
    are discarded and no partial statistics are returned.
    """
    conn = db.connect()
    try:
        row = db.get_task(conn, task_id)
        if row is None:
            raise IndexCheckError(404, "task not found")
        if row["status"] != "succeeded":
            raise IndexCheckError(
                409,
                f"task is {row['status']}; index check requires published "
                "successful results",
            )
    finally:
        conn.close()

    r1_path, r2_path = _resolve_result_paths(task_id, row)

    counts = [0] * len(samples)
    ambiguous = cross_sample = unassigned = total = 0
    try:
        records_r1 = iter_result_records(r1_path, "R1")
        records_r2 = iter_result_records(r2_path, "R2")
        while True:
            rec1 = next(records_r1, None)
            rec2 = next(records_r2, None)
            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                raise FastqFormatError(
                    "R1/R2 record count mismatch: paired results are corrupt"
                )
            total += 1
            (id1, prefix1), (id2, prefix2) = rec1, rec2
            if id1 != id2:
                raise FastqFormatError(
                    f"record {total}: paired read id mismatch R1='{id1}' "
                    f"R2='{id2}'"
                )
            verdict, sample_index = classify_pair(prefix1, prefix2, samples)
            if verdict == "sample":
                counts[sample_index] += 1
            elif verdict == "ambiguous":
                ambiguous += 1
            elif verdict == "cross_sample":
                cross_sample += 1
            else:
                unassigned += 1
    except FastqFormatError:
        raise
    except OSError as exc:
        raise IndexCheckError(500, f"failed reading result file: {exc!r}") from exc

    return {
        "task_id": task_id,
        "total_pairs": total,
        "samples": [
            {"name": spec.name, "count": counts[i]}
            for i, spec in enumerate(samples)
        ],
        "ambiguous": ambiguous,
        "suspected_cross_sample": cross_sample,
        "unassigned": unassigned,
    }
