#!/usr/bin/env bash
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO"
OUT="${VERIFY_OUT:-scratch/verify/truth}"
DB_DIR="${DEMO_DB_DIR:-scratch/demo}"
CONFIRM_SAMPLE="${VERIFY_CONFIRM_SAMPLE:-MO_H2O}"
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

trap 'bin/verify stop > /dev/null 2>&1 || true' EXIT
bin/verify seed --reset > "$OUT/0-seed.json" 2>&1 || { echo "seed failed"; cat "$OUT/0-seed.json"; tail -50 scratch/verify/demo.log 2>/dev/null; exit 1; }
echo "ok      0-seed"

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

if (( ${#failed[@]} )); then
  echo "verify-truth failed: ${failed[*]}"
  exit 1
fi
echo "verify-truth passed"
