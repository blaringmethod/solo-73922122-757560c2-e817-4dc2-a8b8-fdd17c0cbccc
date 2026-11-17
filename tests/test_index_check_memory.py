"""Memory stress test for streaming index checks on very long reads.

Builds a succeeded task whose result records carry sequence/quality lines far
longer than the 256 MiB process cap (default ~100 MiB per line), then runs
run_index_check in a child process and reads its peak RSS from
/proc/<pid>/status (VmHWM). A non-streaming implementation would have to hold
whole lines (or whole files, ~800 MiB here) in memory and blow the cap.

Usage:
    .venv/bin/python -I tests/test_index_check_memory.py [records]
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECORDS = int(sys.argv[1]) if len(sys.argv) > 1 else 10
# One sequence/quality line per record is ~100 MiB; both result files together
# are about 400 MiB of sequence data alone, well above the 256 MiB cap.
LINE_LEN = 100 * 1024 * 1024
RSS_LIMIT_KIB = 120 * 1024

CHILD = textwrap.dedent("""
    import json
    import sys
    sys.path.insert(0, %r)
    from pathlib import Path
    from app import config, db
    from app.indexcheck import SampleSpec, run_index_check

    data_dir = Path(sys.argv[1])
    config.DATA_DIR = data_dir
    config.DB_PATH = data_dir / "tasks.db"
    config.TASKS_DIR = data_dir / "tasks"
    db._initialized = False
    db.init_db()

    samples = [
        SampleSpec("S1", "AAAAAAAA", "CCCCCCCC"),
        SampleSpec("S2", "GGGGGGGG", "TTTTTTTT"),
    ]
    stats = run_index_check("idxmem", samples)
    print("STATS=" + json.dumps(stats, sort_keys=True))
""") % (str(ROOT),)


def write_result_files(base: Path, n: int) -> None:
    out = base / "output"
    out.mkdir(parents=True)
    r1 = out / "R1.fastq"
    r2 = out / "R2.fastq"
    # Two alternating prefixes so both samples get exercised; line bodies are
    # arbitrary padding that the scanner must measure without retaining.
    body_r1 = b"A" * (LINE_LEN - 8)
    body_r2 = b"T" * (LINE_LEN - 8)
    qual_r1 = b"~" * LINE_LEN
    with r1.open("wb") as f1, r2.open("wb") as f2:
        for i in range(n):
            if i % 2 == 0:
                p1 = b"AAAAAAAA"
                p2 = b"CCCCCCCC"
            else:
                p1 = b"GGGGGGGG"
                p2 = b"TTTTTTTT"
            f1.write(b"@long-%010d/1\n" % i + p1 + body_r1 + b"\n+\n"
                     + qual_r1 + b"\n")
            f2.write(b"@long-%010d/1\n" % i + p2 + body_r2 + b"\n+\n"
                     + qual_r1 + b"\n")


def read_vmhwm(pid: int) -> int:
    for line in Path(f"/proc/{pid}/status").read_text().splitlines():
        if line.startswith("VmHWM:"):
            return int(line.split()[1])
    raise RuntimeError("VmHWM not available")


def main() -> int:
    scratch = Path(tempfile.mkdtemp(prefix="fastq-idxmem-",
                                    dir=str(ROOT / "tests" / "output")))
    data_dir = scratch / "data"
    base = data_dir / "tasks" / "idxmem"
    base.mkdir(parents=True)
    print(f"Writing {RECORDS} paired records with {LINE_LEN / (1024**2):.0f} "
          "MiB lines...")
    t0 = time.time()
    write_result_files(base, RECORDS)
    size1 = (base / "output" / "R1.fastq").stat().st_size
    size2 = (base / "output" / "R2.fastq").stat().st_size
    print(f"  inputs ready in {time.time() - t0:.1f}s: "
          f"R1={size1 / (1024**2):.0f} MiB, R2={size2 / (1024**2):.0f} MiB")

    # Register a succeeded task row pointing at the generated results.
    sys.path.insert(0, str(ROOT))
    from app import config, db  # noqa: E402

    config.DATA_DIR = data_dir
    config.DB_PATH = data_dir / "tasks.db"
    config.TASKS_DIR = data_dir / "tasks"
    db._initialized = False
    db.init_db()
    conn = db.connect()
    db.create_task(conn, "idxmem", "R1.fastq", "R2.fastq", size1, size2)
    db.update_task(
        conn, "idxmem", status="succeeded", stage=config.STAGE_DONE,
        progress=1.0, result_r1="output/R1.fastq", result_r2="output/R2.fastq",
    )
    conn.close()

    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    proc = subprocess.Popen(
        [sys.executable, "-I", "-c", CHILD, str(data_dir)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env, text=True,
    )
    t0 = time.time()
    peak = 0
    while proc.poll() is None:
        try:
            peak = max(peak, read_vmhwm(proc.pid))
        except (FileNotFoundError, ProcessLookupError):
            break
        time.sleep(0.1)
    out, err = proc.communicate()
    runtime = time.time() - t0
    print(out.strip())
    if proc.returncode != 0:
        print(err[-4000:])
        print("RESULT: FAIL (worker exited %d)" % proc.returncode)
        return 1

    import json

    stats = json.loads(out.split("STATS=", 1)[1])
    expected = RECORDS // 2
    ok_counts = (
        stats["total_pairs"] == RECORDS
        and stats["samples"][0]["count"] == expected
        and stats["samples"][1]["count"] == RECORDS - expected
        and stats["ambiguous"] == 0
        and stats["suspected_cross_sample"] == 0
        and stats["unassigned"] == 0
    )
    print(f"peak RSS: {peak / 1024:.1f} MiB "
          f"(limit {RSS_LIMIT_KIB / 1024:.0f} MiB), runtime {runtime:.1f}s")
    if not ok_counts:
        print("RESULT: FAIL (unexpected statistics)", stats)
        return 1
    if peak > RSS_LIMIT_KIB:
        print("RESULT: FAIL (memory cap exceeded)")
        return 1
    print("RESULT: PASS")
    shutil.rmtree(scratch, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
