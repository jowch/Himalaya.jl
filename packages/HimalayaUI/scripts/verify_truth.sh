#!/usr/bin/env bash
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO"
OUT="${VERIFY_OUT:-scratch/verify/truth}"
DB_DIR="${DEMO_DB_DIR:-scratch/demo}"
CONFIRM_SAMPLE="${VERIFY_CONFIRM_SAMPLE:-MO_H2O}"
export DEMO_UI_MODE="${DEMO_UI_MODE:-build}"
HELP=packages/HimalayaUI/scripts/verify_truth.mjs
mkdir -p "$OUT"
failed=()

step() {
  local name=$1 ext=$2; shift 2
  if "$@" > "$OUT/$name.$ext" 2>&1; then
    echo "ok      $name"
  else
    echo "FAILED  $name (see $OUT/$name.$ext)"
    head -c 6000 "$OUT/$name.$ext"; echo
    failed+=("$name")
  fi
}

trap 'bin/verify stop > /dev/null 2>&1 || true; bin/verify close > /dev/null 2>&1 || true' EXIT
bin/verify seed --reset --late-trace > "$OUT/0-seed.json" 2>&1 || { echo "seed failed"; cat "$OUT/0-seed.json"; tail -50 scratch/verify/demo.log 2>/dev/null; exit 1; }
echo "ok      0-seed"

step 1-held json bin/verify check-truth --ui --expect fresh
step 1-deliver json bin/verify deliver-late
step 1-rescan log bash -c '
  set -e
  bin/verify open /experiments/1
  for _ in $(seq 1 60); do
    bin/verify api GET /api/experiments/1/problems > "$1/1-problems.json"
    grep -q "no_trace" "$1/1-problems.json" || exit 0
    sleep 1
  done
  echo "opening the experiment never analyzed the delivered traces"; cat "$1/1-problems.json"; exit 1' _ "$OUT"
step 1-fresh json bin/verify check-truth --ui --expect fresh
step 2-triage json bin/verify triage
step 3-triaged json bin/verify check-truth --ui --expect triaged
step 4-confirm json node "$HELP" confirm "$OUT/3-triaged.json" "$CONFIRM_SAMPLE"
step 5-confirmed json bin/verify check-truth --ui --expect triaged
step 5-shows-confirmed json node "$HELP" confirmed "$OUT/5-confirmed.json" "$CONFIRM_SAMPLE"
step 6-reindex log env HIMALAYA_DB_PATH="$(cd "$DB_DIR" && pwd)/himalaya.db" \
  julia --project=packages/HimalayaUI -e '
    using HimalayaUI, DBInterface, Tables
    db = HimalayaUI.open_db()
    id = Int(HimalayaUI._resolve_experiment(db, "Demo (synthetic)").id)
    DBInterface.execute(db, "UPDATE exposures SET trace_hash = NULL, analysis_inputs_hash = NULL WHERE experiment_id = ?", [id])
    HimalayaUI.main(["analyze", "-e", string(id)])
    stale = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM exposures WHERE experiment_id = ? AND status IS NOT ? AND analysis_inputs_hash IS NULL", [id, "rejected"]))).n
    stale == 0 || error("$stale kept exposures were not re-indexed")'
step 7-after-reindex json bin/verify check-truth --ui --expect triaged
step 8-reindex-changed-nothing json node "$HELP" unchanged "$OUT/5-confirmed.json" "$OUT/7-after-reindex.json"
step 10-custom-index log bash -c '
  set -e
  out=$1 help=$2 label=$3
  read -r sid _ < <(node "$help" target "$out/5-confirmed.json" "$label")
  panel="[data-testid=workspace-custom-index]"
  list="[data-testid=workspace-candidates]"
  custom_row() {
    for _ in 1 2 3 4 5; do
      bin/verify snapshot --selector "$list" --raw > "$out/10-$1.yml" || true
      if grep -E "button \".*Custom\"" "$out/10-$1.yml" | grep -q "\[pressed\]"; then [[ $1 == added ]] && return 0
      elif grep -qE "button \".*Custom\"" "$out/10-$1.yml"; then [[ $1 == removed ]] && return 0
      fi
      sleep 1
    done
    echo "the Phase list never showed the custom row $1 (see $out/10-$1.yml)"; return 1
  }
  bin/verify open "/experiments/1/workspace/$sid"
  bin/verify press i
  bin/verify click --role button --name Hexagonal --exact --within "$panel"
  bin/verify click --role button --name "Add Hexagonal" --within "$panel"
  custom_row added
  node "$help" custom "$out/5-confirmed.json" "$label" Hexagonal on
  bin/verify click --role button --name Custom --within "$list"
  custom_row removed
  node "$help" custom "$out/5-confirmed.json" "$label" Hexagonal off' _ "$OUT" "$HELP" "$CONFIRM_SAMPLE"

step 11-peak-clicks log bash -c '
  set -e
  out=$1 help=$2
  twin=$(node "$help" twin "$out/5-confirmed.json") || { echo "$twin"; exit 1; }
  echo "$twin"
  get() { node -e "process.stdout.write(String(JSON.parse(process.argv[1])[process.argv[2]]))" "$twin" "$1"; }
  sid=$(get sample_id) eid=$(get exposure_id) id=$(get id)
  bin/verify open "/experiments/1/workspace/$sid"
  bin/verify click --role button --name "$(get frame)" --exact
  bin/verify click --role button --name "$(get auto)" --exact
  node "$help" twin-state "$eid" "$id" excluded
  bin/verify click --role button --name "$(get excluded)" --exact
  node "$help" twin-state "$eid" "$id" restored
  bin/verify click --role button --name "$(get manual)" --exact
  node "$help" twin-state "$eid" "$id" removed' _ "$OUT" "$HELP"

if (( ${#failed[@]} )); then
  echo "verify-truth failed: ${failed[*]}"
  exit 1
fi
echo "verify-truth passed"
