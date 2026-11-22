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
step 1-late-in-place log bash -c '
  set -e
  out=$1
  sid=$(bin/verify api GET "/api/samples?experiment_id=1" | node -e "let d=\"\";process.stdin.on(\"data\",c=>d+=c).on(\"end\",()=>process.stdout.write(String(JSON.parse(d).body.find(s=>s.name.startsWith(\"MO_lowH2O\")).id)))")
  bin/verify click --selector "[data-testid=workspace-samples] [data-sample-id=\"$sid\"]" --expect-http "404:/api/exposures/[0-9]+/trace"
  plot() { bin/verify count --selector "$1" | node -e "let d=\"\";process.stdin.on(\"data\",c=>d+=c).on(\"end\",()=>{const r=JSON.parse(d);process.stdout.write(r.count+\"|\"+(r.texts[0]??\"\"))})"; }
  before=$(plot "[data-testid=workspace-trace-empty]")
  echo "before delivery: $before"
  [[ $before == "1|Waiting for this image'"'"'s trace." ]] || { echo "the open plot does not say the image is waiting for its trace"; exit 1; }
  bin/verify deliver-late
  bin/verify api POST /api/experiments/1/scan
  for _ in $(seq 1 60); do
    empty=$(plot "[data-testid=workspace-trace-empty]"); markers=$(plot "[data-peak-key]")
    if [[ ${empty%%|*} == 0 && ${markers%%|*} -gt 0 ]]; then echo "after delivery, in place: trace drawn with ${markers%%|*} peak markers"; exit 0; fi
    sleep 1
  done
  echo "after delivery the open workspace still shows: plot message ${empty}, ${markers%%|*} peak markers (no reload)"; exit 1' _ "$OUT"
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
step 7-after-reindex json bin/verify check-truth --ui --reload --expect triaged
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

step 12-bridge-labels log bash -c '
  set -e
  out=$1
  json() { node -e "let d=\"\";process.stdin.on(\"data\",c=>d+=c).on(\"end\",()=>process.stdout.write(String(($1)(JSON.parse(d)))))"; }
  pair=$(bin/verify api GET "/api/samples?experiment_id=1" | json "r=>r.body.slice(0,2).map(s=>s.id).join(\" \")")
  read -r a b <<< "$pair"
  bin/verify api POST "/api/samples/$a/labels" --data "{\"label\":\"37C\"}" > /dev/null
  bin/verify api POST "/api/samples/$a/bridge" --data "{\"to_sample_id\":$b}" > /dev/null
  bin/verify api POST "/api/samples/$b/labels" --data "{\"label\":\"pH 7\"}" > /dev/null
  bin/verify check-truth > "$out/12-bridged.json" || true
  node -e "
    const r = JSON.parse(require(\"fs\").readFileSync(process.argv[1], \"utf8\"));
    const inv = r.invariants.bridged_labels_match;
    if (!inv.checked) { console.log(\"bridged_labels_match did not run: check-truth saw no bridge\"); process.exit(1); }
    if (!inv.ok) { console.log(JSON.stringify(inv.violations, null, 2)); process.exit(1); }
    console.log(\"bridged_labels_match checked and held\");" "$out/12-bridged.json"
  for s in $a $b; do
    got=$(bin/verify api GET "/api/experiments/1/labels" | json "r=>r.body.filter(l=>l.sample_id===$s&&l.state===\"confirmed\").map(l=>l.label).sort().join(\",\")")
    [[ $got == *37C* && $got == *"pH 7"* ]] || { echo "sample $s has labels [$got], not both 37C and pH 7"; exit 1; }
  done
  echo "both bridged samples carry 37C and pH 7"' _ "$OUT"

step 13-bridge-call json node "$HELP" bridge-call
step 13-calls-match log bash -c '
  out=$1
  bin/verify check-truth > "$out/13-bridged.json" || true
  node -e "
    const r = JSON.parse(require(\"fs\").readFileSync(process.argv[1], \"utf8\"));
    const inv = r.invariants.bridged_calls_match;
    if (!inv.checked) { console.log(\"bridged_calls_match did not run: check-truth saw no bridge\"); process.exit(1); }
    if (!inv.ok) { console.log(JSON.stringify(inv.violations, null, 2)); process.exit(1); }
    console.log(\"bridged_calls_match checked and held\");" "$out/13-bridged.json"' _ "$OUT"
if (( ${#failed[@]} )); then
  echo "verify-truth failed: ${failed[*]}"
  exit 1
fi
echo "verify-truth passed"
