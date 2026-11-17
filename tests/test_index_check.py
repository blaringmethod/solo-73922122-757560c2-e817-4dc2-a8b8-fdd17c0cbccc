"""Functional self-tests for the index-check (索引核查) endpoint.

Run with the project virtualenv:
    .venv/bin/python -I tests/test_index_check.py
Each run gets an isolated FASTQ_DATA_DIR under tests/output.
"""

from __future__ import annotations

import asyncio
import hashlib
import io
import json
import os
import shutil
import sys
import tempfile
import threading
import time
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import httpx  # noqa: E402

from app import config, db, worker  # noqa: E402
from app.main import app  # noqa: E402

RESULTS: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))


def make_record(rid: str, seq: bytes, quality: bytes | None = None,
                suffix: str = "") -> bytes:
    if quality is None:
        quality = b"~" * len(seq)
    return b"\n".join((b"@" + rid.encode() + suffix.encode(), seq, b"+", quality)) + b"\n"


def fastq(records: list[bytes]) -> io.BytesIO:
    return io.BytesIO(b"".join(records))


async def wait_terminal(client: httpx.AsyncClient, task_id: str,
                        timeout: float = 30.0) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        body = (await client.get(f"/tasks/{task_id}")).json()
        if body["status"] in ("succeeded", "failed", "cancelled"):
            return body
        await asyncio.sleep(0.05)
    raise AssertionError(f"task {task_id} did not finish in {timeout}s")


async def upload(client: httpx.AsyncClient, r1: io.BytesIO, r2: io.BytesIO) -> str:
    r1.seek(0)
    r2.seek(0)
    r = await client.post(
        "/tasks",
        files={"r1": ("r1.fastq", r1, "application/octet-stream"),
               "r2": ("r2.fastq", r2, "application/octet-stream")},
    )
    assert r.status_code == 202, r.text
    return r.json()["id"]


SAMPLES = [
    {"name": "S1", "r1_index": "AAAAAAAA", "r2_index": "CCCCCCCC"},
    {"name": "S2", "r1_index": "GGGGGGGG", "r2_index": "TTTTTTTT"},
    # S3 is 2 mismatches away from S1 on each end, so 1mm balls can overlap
    # (for the tie case) while unique pointing stays possible (for crossover).
    {"name": "S3", "r1_index": "AAAAAATT", "r2_index": "CCCCCCTA"},
]

# (r1 prefix, r2 prefix, expected verdict) — pairs 1..9
PAIRS = [
    (b"AAAAAAAA", b"CCCCCCCC", "S1"),          # 1 exact S1
    (b"AAAAAAAT", b"CCCCCTGG", "unassigned"),  # 2 R1 S1/S3, R2: 2mm from S1/S3
    (b"GGGGGGGG", b"CCCCCCCG", "cross"),       # 3 R1 unique->S2, R2 unique->S1
    (b"GGGGGGGG", b"TTTTTTTT", "S2"),          # 4 exact S2
    (b"GGGGGGGT", b"TTTTTTTG", "S2"),          # 5 1mm+1mm -> S2 total 2
    (b"AAAAAAAT", b"CCCCCCTC", "ambiguous"),   # 6 tie at total 2
    (b"AAAAAAAN", b"CCCCCCCC", "S1"),          # 7 non-ACGT fallback, 1mm -> S1
    (b"AAAAAAA", b"CCCCCCC", "unassigned"),    # 8 7-base prefixes match nothing
    (b"NNNNNNNN", b"NNNNNNNN", "unassigned"),  # 9 neither end has a candidate
]
# Pair 2: R1 has two candidates (S1, S3, each 1mm) but R2 matches neither ->
# not a unique/unique-different pair, so it must be 'unassigned', not crossover.
# Pair 6: observed R1 is 1mm from both S1 and S3 R1; observed R2 is 1mm from
# both R2 indices -> both legal combos total 2 mismatches -> tie -> ambiguous.


def build_pair_files() -> tuple[io.BytesIO, io.BytesIO]:
    r1_records, r2_records = [], []
    for i, (p1, p2, _verdict) in enumerate(PAIRS, start=1):
        # Only pad when the prefix is already a full index, so sub-8-base
        # prefixes stay short and must match no sample.
        tail1 = b"ACGTACGT" * 8 if len(p1) >= 8 else b""
        tail2 = b"TGCATGCA" * 8 if len(p2) >= 8 else b""
        r1_records.append(make_record(f"read-{i}", p1 + tail1))
        r2_records.append(make_record(f"read-{i}", p2 + tail2))
    return fastq(r1_records), fastq(r2_records)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


async def main_async() -> int:
    tmp_root = Path(tempfile.mkdtemp(prefix="fastq-idx-",
                                     dir=str(ROOT / "tests" / "output")))
    data_dir = tmp_root / "data"
    os.environ["FASTQ_DATA_DIR"] = str(data_dir)
    config.DATA_DIR = data_dir
    config.DB_PATH = data_dir / "tasks.db"
    config.TASKS_DIR = data_dir / "tasks"
    config.PROGRESS_UPDATE_EVERY_BYTES = 64 * 1024

    db._initialized = False
    db.init_db()
    worker.recover_on_startup()

    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(transport=transport, base_url="http://test",
                                 timeout=60) as client:
        # ---- 1. happy path: full classification matrix
        tid = await upload(client, *build_pair_files())
        body = await wait_terminal(client, tid)
        record("fixture task succeeds", body["status"] == "succeeded",
               body.get("message", ""))

        r = await client.post(f"/tasks/{tid}/index-check", json={"samples": SAMPLES})
        record("index check returns 200", r.status_code == 200, str(r.status_code))
        stats = r.json()
        by_name = {item["name"]: item["count"] for item in stats["samples"]}
        record("total pair count", stats["total_pairs"] == len(PAIRS),
               str(stats["total_pairs"]))
        record("S1 count (exact + non-ACGT fallback)", by_name.get("S1") == 2,
               str(by_name))
        record("S2 count (exact + dual 1mm)", by_name.get("S2") == 2, str(by_name))
        record("S3 count stays zero", by_name.get("S3") == 0, str(by_name))
        record("ambiguous count", stats["ambiguous"] == 1, str(stats["ambiguous"]))
        record("suspected crossover count",
               stats["suspected_cross_sample"] == 1,
               str(stats["suspected_cross_sample"]))
        record("unassigned count", stats["unassigned"] == 3, str(stats["unassigned"]))
        accounted = (sum(by_name.values()) + stats["ambiguous"]
                     + stats["suspected_cross_sample"] + stats["unassigned"])
        record("all pairs accounted exactly once", accounted == stats["total_pairs"],
               str(accounted))
        record("all sheet samples echoed in order",
               [item["name"] for item in stats["samples"]] == ["S1", "S2", "S3"])

        # ---- 2. input validation (pydantic -> 422, no file access)
        async def expect_422(label: str, payload: dict) -> None:
            rr = await client.post(f"/tasks/{tid}/index-check", json=payload)
            record(label, rr.status_code == 422, str(rr.status_code))

        one = SAMPLES[:1]
        await expect_422("fewer than 2 samples rejected", {"samples": one})
        many = [
            {"name": f"X{i}", "r1_index": "AAAAAAAA", "r2_index": "CCCCCCCC"}
            for i in range(33)
        ]
        await expect_422("more than 32 samples rejected", {"samples": many})
        bad_char = [
            {"name": "A", "r1_index": "AAAAAAAN", "r2_index": "CCCCCCCC"},
            {"name": "B", "r1_index": "GGGGGGGG", "r2_index": "TTTTTTTT"},
        ]
        await expect_422("non-ACGT index rejected", {"samples": bad_char})
        bad_len = [
            {"name": "A", "r1_index": "AAAAAAA", "r2_index": "CCCCCCCC"},
            {"name": "B", "r1_index": "GGGGGGGG", "r2_index": "TTTTTTTT"},
        ]
        await expect_422("7-char index rejected", {"samples": bad_len})
        dup = [
            {"name": "DUP", "r1_index": "AAAAAAAA", "r2_index": "CCCCCCCC"},
            {"name": "DUP", "r1_index": "GGGGGGGG", "r2_index": "TTTTTTTT"},
        ]
        await expect_422("duplicate sample names rejected", {"samples": dup})
        blank = [
            {"name": "  ", "r1_index": "AAAAAAAA", "r2_index": "CCCCCCCC"},
            {"name": "B", "r1_index": "GGGGGGGG", "r2_index": "TTTTTTTT"},
        ]
        await expect_422("blank sample name rejected", {"samples": blank})
        await expect_422("empty samples list rejected", {"samples": []})

        # ---- 3. task state guards
        unknown = await client.post(
            f"/tasks/{'0' * 32}/index-check", json={"samples": SAMPLES})
        record("unknown task -> 404", unknown.status_code == 404,
               str(unknown.status_code))

        bad = await upload(
            client,
            fastq([make_record("x", b"AAAA")]),
            fastq([make_record("x", b"AAAA", b"~~")]),  # seq/qual length mismatch
        )
        failed = await wait_terminal(client, bad)
        record("malformed fixture fails", failed["status"] == "failed",
               failed.get("message", ""))
        rr = await client.post(f"/tasks/{bad}/index-check",
                               json={"samples": SAMPLES})
        record("index check on failed task -> 409", rr.status_code == 409,
               str(rr.status_code))

        # ---- 3b. publish-window race: 'succeeded' before files land
        race_id = "race" + uuid.uuid4().hex[:8]
        race_out = config.TASKS_DIR / race_id / "output"
        race_out.mkdir(parents=True)
        conn = db.connect()
        db.create_task(conn, race_id, "R1.fastq", "R2.fastq", 1, 1)
        db.update_task(conn, race_id, status="succeeded",
                       stage=config.STAGE_DONE, progress=1.0,
                       result_r1="output/R1.fastq", result_r2="output/R2.fastq")
        conn.close()
        payload_r1 = make_record("r1x", b"AAAAAAAA" + b"GGGG",
                                 b"~" * 12)
        payload_r2 = make_record("r1x", b"CCCCCCCC" + b"AAAA",
                                 b"!" * 12)

        def publish_later() -> None:
            time.sleep(0.3)  # after the check has started waiting
            (race_out / "R1.fastq").write_bytes(payload_r1)
            (race_out / "R2.fastq").write_bytes(payload_r2)

        threading.Thread(target=publish_later, daemon=True).start()
        rr = await client.post(f"/tasks/{race_id}/index-check",
                               json={"samples": SAMPLES[:2]})
        record("check waits out the publish rename window",
               rr.status_code == 200
               and rr.json()["samples"][0]["count"] == 1, rr.text[:120])

        # ---- 4. result files are never modified by a check
        out_r1 = config.TASKS_DIR / tid / "output" / "R1.fastq"
        out_r2 = config.TASKS_DIR / tid / "output" / "R2.fastq"
        before = (sha256(out_r1), sha256(out_r2))
        await client.post(f"/tasks/{tid}/index-check", json={"samples": SAMPLES})
        after = (sha256(out_r1), sha256(out_r2))
        record("successful check leaves result bytes untouched",
               before == after)

        # ---- 5. missing / corrupt results -> whole error, no partial stats
        missing_dir = config.TASKS_DIR / tid / "output"
        original_r1 = out_r1.read_bytes()
        original_r2 = out_r2.read_bytes()
        try:
            out_r1.unlink()
            rr = await client.post(f"/tasks/{tid}/index-check",
                                   json={"samples": SAMPLES})
            record("missing result -> 422 whole error", rr.status_code == 422,
                   str(rr.status_code))
            record("missing result gives no partial stats",
                   "total_pairs" not in rr.text)
        finally:
            out_r1.write_bytes(original_r1)

        def restore() -> None:
            out_r1.write_bytes(original_r1)
            out_r2.write_bytes(original_r2)

        # truncated final record
        try:
            out_r1.write_bytes(original_r1[:-40])
            rr = await client.post(f"/tasks/{tid}/index-check",
                                   json={"samples": SAMPLES})
            record("truncated record -> 422", rr.status_code == 422, rr.text[:120])
            record("truncation error is whole-request",
                   "total_pairs" not in rr.text)
        finally:
            restore()

        # seq/quality length mismatch introduced into stored result
        try:
            # rewrite the last record's quality line to a shorter value
            lines = original_r1.splitlines(keepends=True)
            lines[-1] = b"~~\n"
            out_r1.write_bytes(b"".join(lines))
            rr = await client.post(f"/tasks/{tid}/index-check",
                                   json={"samples": SAMPLES})
            record("corrupt seq/qual length -> 422", rr.status_code == 422,
                   rr.text[:120])
        finally:
            restore()

        # paired read id mismatch between R1 and R2
        try:
            lines = original_r1.splitlines(keepends=True)
            lines[0] = b"@different-id\n"
            out_r1.write_bytes(b"".join(lines))
            rr = await client.post(f"/tasks/{tid}/index-check",
                                   json={"samples": SAMPLES})
            record("paired id mismatch -> 422", rr.status_code == 422, rr.text[:120])
        finally:
            restore()

        # record count mismatch: extra record only on R2
        try:
            extra = make_record("read-extra", b"AAAAAAAA" + b"ACGT",
                                b"~" * 12)
            out_r2.write_bytes(original_r2 + extra)
            rr = await client.post(f"/tasks/{tid}/index-check",
                                   json={"samples": SAMPLES})
            record("R1/R2 count mismatch -> 422", rr.status_code == 422,
                   rr.text[:120])
        finally:
            restore()

        # failed checks must not have mutated the files
        record("failed checks leave result bytes untouched",
               (sha256(out_r1), sha256(out_r2)) == before)

        # ---- 6. long reads: only the 8-base prefix is retained in memory
        long_seq = b"AAAAAAAA" + b"G" * (2 * 1024 * 1024)
        long_r1 = fastq([
            make_record("L1", long_seq),
            make_record("L2", b"GGGGGGGG" + b"T" * (2 * 1024 * 1024)),
        ])
        long_r2 = fastq([
            make_record("L1", b"CCCCCCCC" + b"A" * (2 * 1024 * 1024)),
            make_record("L2", b"TTTTTTTT" + b"C" * (2 * 1024 * 1024)),
        ])
        lid = await upload(client, long_r1, long_r2)
        body = await wait_terminal(client, lid)
        record("long-read task succeeds", body["status"] == "succeeded",
               body.get("message", ""))
        rr = await client.post(f"/tasks/{lid}/index-check",
                               json={"samples": SAMPLES[:2]})
        record("long-read index check 200", rr.status_code == 200, rr.text[:120])
        long_stats = rr.json()
        long_by_name = {item["name"]: item["count"]
                        for item in long_stats["samples"]}
        record("long reads classified from prefixes only",
               long_stats["total_pairs"] == 2
               and long_by_name == {"S1": 1, "S2": 1},
               json.dumps(long_stats))

        # CRLF line endings and no trailing newline on the final line
        crlf_r1 = io.BytesIO(
            b"@z/1\r\nAAAAAAAACCCC\r\n+\r\n" + b"~" * 12 + b"\r\n")
        crlf_r2 = io.BytesIO(
            b"@z/2\r\nCCCCCCCCGGGG\r\n+\r\n" + b"!" * 12 + b"\r\n")
        cid = await upload(client, crlf_r1, crlf_r2)
        body = await wait_terminal(client, cid)
        record("CRLF fixture succeeds", body["status"] == "succeeded",
               body.get("message", ""))
        rr = await client.post(f"/tasks/{cid}/index-check",
                               json={"samples": SAMPLES[:2]})
        record("CRLF prefixes parsed correctly", rr.status_code == 200
               and rr.json()["samples"][0]["count"] == 1, rr.text[:120])

    shutil.rmtree(tmp_root, ignore_errors=True)
    failed = [name for name, ok, _ in RESULTS if not ok]
    print("\n%d/%d checks passed" % (len(RESULTS) - len(failed), len(RESULTS)))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(asyncio.run(main_async()))
