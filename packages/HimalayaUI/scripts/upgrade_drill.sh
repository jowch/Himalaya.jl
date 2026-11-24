#!/usr/bin/env bash
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO"
BASE_REF="${UPGRADE_BASE:-origin/main}"
ROOT="${UPGRADE_DIR:-scratch/upgrade}"
PORT="${UPGRADE_BASE_PORT:-8197}"
mkdir -p "$ROOT"
ROOT="$(cd "$ROOT" && pwd)"
WT="$ROOT/base"
DEMO="$ROOT/demo-experiment"
DB_DIR="$ROOT/db"
OUT="$ROOT/out"
JL="julia --threads=4,1"
HISTORY=packages/HimalayaUI/scripts/upgrade_history.mjs
HELP_TRUTH=packages/HimalayaUI/scripts/verify_truth.mjs

base_pid=""
cleanup() {
  [[ -n $base_pid ]] && kill "$base_pid" 2>/dev/null || true
  DEMO_DIR="$DEMO" DEMO_DB_DIR="$DB_DIR" bin/verify stop > /dev/null 2>&1 || true
  bin/verify close > /dev/null 2>&1 || true
  git worktree remove --force "$WT" > /dev/null 2>&1 || true
}
trap cleanup EXIT

rm -rf "$DB_DIR" "$OUT" "$DEMO"
git worktree remove --force "$WT" > /dev/null 2>&1 || rm -rf "$WT"
mkdir -p "$DB_DIR" "$OUT"
sha=$(git rev-parse --verify "$BASE_REF^{commit}")
echo "base: $BASE_REF ($sha)"
git worktree add --detach "$WT" "$sha" > /dev/null

$JL --project=packages/HimalayaUI packages/HimalayaUI/scripts/make_demo_data.jl "$DEMO" > "$OUT/demo-data.log" 2>&1
$JL --project="$WT/packages/HimalayaUI" -e "using Pkg; Pkg.develop(path=\"$WT\"); Pkg.instantiate()" > "$OUT/base-install.log" 2>&1

(cd "$WT" && HIMALAYA_DB_PATH="$DB_DIR/himalaya.db" exec $JL --project=packages/HimalayaUI \
  -e 'using HimalayaUI; main(ARGS)' -- serve --port "$PORT") > "$OUT/base-backend.log" 2>&1 &
base_pid=$!
for _ in $(seq 1 300); do
  curl -fsS "http://127.0.0.1:$PORT/api/experiments" > /dev/null 2>&1 && break
  kill -0 "$base_pid" 2>/dev/null || { echo "the base backend exited; see $OUT/base-backend.log"; exit 1; }
  sleep 1
done
(cd "$WT" && $JL --project=packages/HimalayaUI packages/HimalayaUI/scripts/demo_ingest.jl "$DEMO" --port "$PORT") > "$OUT/base-ingest.log" 2>&1
review="$WT/$HISTORY"
[[ -f $review ]] || review="$HISTORY"
echo "review script: ${review#$REPO/}"
node "$review" review "$PORT" "$OUT/calls-before.json"
kill "$base_pid"; wait "$base_pid" 2>/dev/null || true; base_pid=""

export DEMO_DIR="$DEMO" DEMO_DB_DIR="$DB_DIR" DEMO_UI_MODE="${DEMO_UI_MODE:-build}"
bin/verify seed > "$OUT/seed.json" 2>&1 || { echo "the head app did not come up on the upgraded DB"; cat "$OUT/seed.json"; tail -50 "$DB_DIR/backend.log" 2>/dev/null; exit 1; }
failed=0
node "$HISTORY" check "${DEMO_API_PORT:-8080}" "$OUT/calls-before.json" > "$OUT/calls-after.json" || { echo "FAILED  calls changed (see $OUT/calls-after.json)"; cat "$OUT/calls-after.json"; failed=1; }
bin/verify check-truth --ui --gate invariants > "$OUT/check-truth.json" || { node "$HELP_TRUTH" broken "$OUT/check-truth.json"; echo "FAILED  check-truth invariants on the upgraded DB (see $OUT/check-truth.json)"; failed=1; }
node -e '
  const r = JSON.parse(require("fs").readFileSync(process.argv[1], "utf8"));
  const inv = r.invariants ?? {};
  const unchecked = Object.entries(inv).filter(([k, v]) => !v.checked && !["bridged_labels_match", "bridged_calls_match"].includes(k)).map(([k]) => k);
  const broken = Object.entries(inv).filter(([, v]) => v.ok === false).map(([k, v]) => ({ invariant: k, violations: v.violations }));
  const stale = (r.samples ?? []).filter((s) => !s.ok).map((s) => ({ sample: s.sample, problems: s.problems }));
  if (stale.length) console.log("auto calls that differ from TRUTH (analyses the base wrote; not an upgrade failure):\n" + JSON.stringify(stale, null, 2));
  const bad = broken.length || unchecked.length || !r.traces?.ok || (r.strays ?? []).length;
  if (bad) { console.log(JSON.stringify({ broken, unchecked, traces: r.traces, strays: r.strays }, null, 2)); process.exit(1); }
  console.log("every check-truth invariant was checked and held on the upgraded DB");' "$OUT/check-truth.json" \
  || { echo "FAILED  check-truth invariants on the upgraded DB (see $OUT/check-truth.json)"; failed=1; }
missing=$(node -e 'process.stdout.write(JSON.parse(require("fs").readFileSync(process.argv[1],"utf8")).reviewed[0].sample)' "$OUT/calls-before.json")
sid=$(bin/verify api GET "/api/samples?experiment_id=1" | node -e 'let d="";process.stdin.on("data",c=>d+=c).on("end",()=>process.stdout.write(String(JSON.parse(d).body.find(s=>s.name===process.argv[1])?.id??"")))' "$missing")
bin/verify open "/experiments/1/workspace/$sid" > "$OUT/open-not-found.json"
marked=$(bin/verify count --selector "[data-testid=workspace-candidates] [data-unsupported][aria-pressed=true]" | node -e 'let d="";process.stdin.on("data",c=>d+=c).on("end",()=>process.stdout.write(String(JSON.parse(d).count)))')
if [[ $marked -ge 1 ]]; then echo "$missing shows its call ticked and marked Not found"
else echo "FAILED  $missing's Assignment list has no ticked Not found row, so rail_matches_call never saw one"; failed=1; fi
(( failed )) && exit 1
echo "upgrade drill passed: $BASE_REF's DB with reviews opens on this head with the same calls and every check-truth invariant"
