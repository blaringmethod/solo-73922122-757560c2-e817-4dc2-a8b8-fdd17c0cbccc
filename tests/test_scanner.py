"""Unit tests for the bounded chunk-streaming FASTQ scanner.

Covers line-ending and truncation corner cases that are cheap to construct
but easy to get wrong in a hand-rolled chunk reader.
"""

from __future__ import annotations

import io
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from app.fastq import FastqFormatError  # noqa: E402
from app.indexcheck import _ChunkReader, iter_result_records  # noqa: E402

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    if not ok:
        FAILURES.append(name)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))


def reader(data: bytes, chunk: int = 4) -> _ChunkReader:
    return _ChunkReader(io.BytesIO(data), chunk_size=chunk)


def main() -> int:
    # 1. plain LF, content length excludes the newline
    r = reader(b"abc\n")
    line, size, eof = r.read_line(cap=16)
    check("basic LF line", line == b"abc" and size == 3 and not eof)
    line, size, eof = r.read_line(cap=16)
    check("clean EOF", line == b"" and size == 0 and eof)

    # 2. CRLF counted once, including when split across chunk boundaries
    for chunk in (1, 2, 3, 4, 5, 7, 64):
        r = reader(b"ab\r\ncd\r\n", chunk=chunk)
        l1, s1, e1 = r.read_line(cap=16)
        l2, s2, e2 = r.read_line(cap=16)
        check(f"CRLF split across chunks (chunk={chunk})",
              (l1, s1) == (b"ab", 2) and (l2, s2) == (b"cd", 2) and not e2)

    # 3. lone CR at EOF is content (no LF follows), CR mid-line is content
    r = reader(b"a\rb\r\n")
    l1, s1, _ = r.read_line(cap=16)
    check("mid-line CR retained, CRLF CR stripped", l1 == b"a\rb" and s1 == 3)

    r = reader(b"ab\r", chunk=2)
    line, size, eof = r.read_line(cap=16)
    check("lone CR at EOF retained as content",
          line == b"ab\r" and size == 3 and eof)

    # 4. final line without trailing newline still parses
    r = reader(b"abc")
    line, size, eof = r.read_line(cap=16)
    check("EOF-terminated line", line == b"abc" and size == 3 and eof)

    # 5. long line: only prefix retained, full length still measured
    big = b"X" * 100 + b"\n"
    r = reader(big, chunk=8)
    prefix, size, eof = r.read_line(keep_prefix=8)
    check("long line measured, prefix kept",
          prefix == b"X" * 8 and size == 100 and not eof)

    # 6. cap mode enforces the limit and retains full content
    r = reader(b"123456789\n", chunk=3)
    try:
        r.read_line(cap=4)
        check("line over cap rejected", False)
    except FastqFormatError:
        check("line over cap rejected", True)
    r = reader(b"1234\n", chunk=3)
    line, size, _ = r.read_line(cap=4)
    check("line at cap accepted", line == b"1234" and size == 4)

    # 7. keep_prefix with CRLF and chunk splits
    r = reader(b"ABCDEFGHIJ\r\n", chunk=3)
    prefix, size, _ = r.read_line(keep_prefix=8)
    check("prefix mode strips CRLF", prefix == b"ABCDEFGH" and size == 10)

    # 8. full record iterator: no trailing newline on quality line is valid
    data = b"@r/1\nACGTACGT\n+\n~~~~~~~~"  # final line, no LF
    tmp = ROOT / "tests" / "output" / "scanner_tmp.fastq"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    tmp.write_bytes(data)
    recs = list(iter_result_records(tmp, "R1"))
    check("record without trailing newline",
          len(recs) == 1 and recs[0] == ("r", b"ACGTACGT"))

    # 9. truncation variants all raise
    for label, bad in (
        ("missing plus", b"@r\nACGT\n"),
        ("missing quality", b"@r\nACGT\n+\n"),
        ("bad header", b"r\nACGT\n+\n~~~~\n"),
        ("bad plus", b"@r\nACGT\n?\n~~~~\n"),
        ("length mismatch", b"@r\nACGT\n+\n~~\n"),
        ("empty file", b""),
    ):
        tmp.write_bytes(bad)
        try:
            list(iter_result_records(tmp, "R1"))
            check(f"{label} rejected", False)
        except FastqFormatError:
            check(f"{label} rejected", True)

    tmp.unlink(missing_ok=True)
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
