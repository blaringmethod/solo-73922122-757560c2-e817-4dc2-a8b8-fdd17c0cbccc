#!/usr/bin/env bash
# End-to-end demo against a running service (default http://localhost:8000).
# Shows: upload -> progress polling -> download verification, plus cancellation
# and a rejected (malformed) upload.
set -euo pipefail

BASE_URL="${BASE_URL:-http://localhost:8000}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

need() { command -v "$1" >/dev/null || { echo "missing tool: $1" >&2; exit 1; }; }
need curl
need python3

echo "== health =="
curl -fsS "$BASE_URL/health"
echo

mk() { # file id seq side
  printf '@%s/%s desc\n%s\n+\n%s\n' "$2" "$4" "$3" "$(printf '~%.0s' $(seq ${#3}))" >> "$1"
}

echo "== generating shuffled R1/R2 =="
R1="$WORK/R1.fastq"; R2="$WORK/R2.fastq"
for id in alpha beta gamma delta epsilon; do
  mk "$R1" "$id" "ACGTACGTAC" 1
done
for id in epsilon gamma alpha delta beta; do   # deliberately shuffled
  mk "$R2" "$id" "TGCATGCATG" 2
done

echo "== upload =="
TASK=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$R1;filename=R1.fastq" \
  -F "r2=@$R2;filename=R2.fastq")
echo "$TASK"
ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$TASK")

echo "== poll until terminal =="
while :; do
  BODY=$(curl -fsS "$BASE_URL/tasks/$ID")
  STATUS=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["status"])' <<<"$BODY")
  echo "$BODY"
  case "$STATUS" in
    succeeded|failed|cancelled) break ;;
  esac
  sleep 0.5
done

echo "== download and verify order =="
curl -fsS "$BASE_URL/tasks/$ID/download/r1" -o "$WORK/out_R1.fastq"
curl -fsS "$BASE_URL/tasks/$ID/download/r2" -o "$WORK/out_R2.fastq"
python3 - "$WORK/out_R1.fastq" "$WORK/out_R2.fastq" <<'PY'
import sys
ids = []
for path in sys.argv[1:]:
    got = []
    with open(path, "rb") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            got.append(header[1:].split()[0].decode().rsplit("/", 1)[0])
            fh.readline(); fh.readline(); fh.readline()
    ids.append(got)
print("R1 ids:", ids[0])
print("R2 ids:", ids[1])
assert ids[0] == ids[1] == ["alpha", "beta", "gamma", "delta", "epsilon"]
print("OK: paired in R1 order")
PY

echo "== malformed input is rejected as a whole task =="
printf '@x\nACGT\n+\n~~\n' > "$WORK/bad_R2.fastq"     # seq 4 != qual 2
BAD=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$R1;filename=R1.fastq" \
  -F "r2=@$WORK/bad_R2.fastq;filename=R2.fastq")
BAD_ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$BAD")
sleep 1
curl -fsS "$BASE_URL/tasks/$BAD_ID"
echo
echo "(download attempt on failed task -> HTTP 409)"
curl -s -o /dev/null -w '%{http_code}\n' "$BASE_URL/tasks/$BAD_ID/download/r1"

echo "== cancellation returns immediately, task never publishes =="
BIG_R1="$WORK/big_R1.fastq"; BIG_R2="$WORK/big_R2.fastq"
python3 - "$BIG_R1" "$BIG_R2" <<'PY'
import sys
n = 200_000
with open(sys.argv[1], "wb") as f1, open(sys.argv[2], "wb") as f2:
    for i in range(n):
        f1.write(b"@r%09d/1\n" % i + b"A" * 40 + b"\n+\n" + b"~" * 40 + b"\n")
    for i in range(n - 1, -1, -1):
        f2.write(b"@r%09d/2\n" % i + b"T" * 40 + b"\n+\n" + b"!" * 40 + b"\n")
PY
CANCEL=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$BIG_R1;filename=big_R1.fastq" \
  -F "r2=@$BIG_R2;filename=big_R2.fastq")
CANCEL_ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$CANCEL")
curl -fsS -X POST "$BASE_URL/tasks/$CANCEL_ID/cancel"
echo
sleep 2
curl -fsS "$BASE_URL/tasks/$CANCEL_ID"
echo

echo "== cancelling the already succeeded task -> HTTP 409 =="
curl -s -o /dev/null -w '%{http_code}\n' -X POST "$BASE_URL/tasks/$ID/cancel"

echo "== index check (索引核查): controlled index prefixes =="
# Two samples on the sheet; each read carries an 8-base ACGT index prefix.
IC_R1="$WORK/ic_R1.fastq"; IC_R2="$WORK/ic_R2.fastq"
python3 - "$IC_R1" "$IC_R2" <<'PY'
import sys
r1, r2 = sys.argv[1], sys.argv[2]
pairs = [
    # id, R1 prefix+body, R2 prefix+body, expected classification
    ("p1", b"AAAAAAAA" + b"G" * 32, b"CCCCCCCC" + b"A" * 32),  # exact sample-1
    ("p2", b"AAAAAAAG" + b"G" * 32, b"CCCCCCCA" + b"A" * 32),  # 1mm+1mm sample-1
    ("p3", b"GGGGGGGG" + b"T" * 32, b"CCCCCCCC" + b"A" * 32),  # crossover
    ("p4", b"NNNNNNNN" + b"C" * 32, b"NNNNNNNN" + b"G" * 32),  # unassigned
]
with open(r1, "wb") as f1, open(r2, "wb") as f2:
    for rid, s1, s2 in pairs:
        f1.write(b"@%s/1\n%s\n+\n%s\n" % (rid.encode(), s1, b"~" * len(s1)))
        f2.write(b"@%s/2\n%s\n+\n%s\n" % (rid.encode(), s2, b"!" * len(s2)))
PY
IC=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$IC_R1;filename=ic_R1.fastq" \
  -F "r2=@$IC_R2;filename=ic_R2.fastq")
IC_ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$IC")
while :; do
  BODY=$(curl -fsS "$BASE_URL/tasks/$IC_ID")
  STATUS=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["status"])' <<<"$BODY")
  case "$STATUS" in
    succeeded|failed|cancelled) break ;;
  esac
  sleep 0.3
done
echo "index-check task status: $STATUS"
curl -fsS -X POST "$BASE_URL/tasks/$IC_ID/index-check" \
  -H 'Content-Type: application/json' \
  -d '{"samples":[
        {"name":"sample-1","r1_index":"AAAAAAAA","r2_index":"CCCCCCCC"},
        {"name":"sample-2","r1_index":"GGGGGGGG","r2_index":"TTTTTTTT"}
      ]}'
echo

echo "== index check input guards: bad sheet -> 422, failed task -> 409 =="
curl -s -o /dev/null -w 'bad index chars -> %{http_code}\n' -X POST \
  "$BASE_URL/tasks/$IC_ID/index-check" -H 'Content-Type: application/json' \
  -d '{"samples":[{"name":"a","r1_index":"AAAAAAAN","r2_index":"CCCCCCCC"},
                  {"name":"b","r1_index":"GGGGGGGG","r2_index":"TTTTTTTT"}]}'
curl -s -o /dev/null -w 'check on failed task -> %{http_code}\n' -X POST \
  "$BASE_URL/tasks/$BAD_ID/index-check" -H 'Content-Type: application/json' \
  -d '{"samples":[{"name":"a","r1_index":"AAAAAAAA","r2_index":"CCCCCCCC"},
                  {"name":"b","r1_index":"GGGGGGGG","r2_index":"TTTTTTTT"}]}'

echo "ALL DEMO STEPS DONE"
