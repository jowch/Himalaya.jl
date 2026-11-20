#!/usr/bin/env bash
# Live check of the synthetic demo against its TRUTH.csv, through the real app:
# a fresh demo (regenerated, ingested through the New-experiment funnel), after
# triage drops TRUTH's bad frames, and after `himalaya analyze` re-runs the
# pipeline, which must not change what any view shows. `make verify-truth` and CI run this. Each step's JSON lands in
# $VERIFY_OUT; the script exits non-zero if any step fails.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO"
OUT="${VERIFY_OUT:-scratch/verify/truth}"
DB_DIR="${DEMO_DB_DIR:-scratch/demo}"
mkdir -p "$OUT"
failed=()

step() {
  local name=$1; shift
  if "$@" > "$OUT/$name.json" 2>&1; then
    echo "ok      $name"
  else
    echo "FAILED  $name (see $OUT/$name.json)"
    head -c 6000 "$OUT/$name.json"; echo
    failed+=("$name")
  fi
}

bin/verify seed --reset > "$OUT/0-seed.json" 2>&1 || { echo "seed failed"; cat "$OUT/0-seed.json"; tail -50 scratch/verify/demo.log 2>/dev/null; exit 1; }
echo "ok      0-seed"
trap 'bin/verify stop > /dev/null 2>&1 || true' EXIT

step 1-fresh bin/verify check-truth --ui
step 2-triage bin/verify triage
step 3-triaged bin/verify check-truth --ui
step 4-analyze env HIMALAYA_DB_PATH="$(cd "$DB_DIR" && pwd)/himalaya.db" \
  julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- analyze -e "Demo (synthetic)"
step 5-after-analyze bin/verify check-truth --ui
step 6-analyze-changed-nothing node -e '
  const [a, b] = process.argv.slice(1).map((f) => JSON.parse(require("fs").readFileSync(f, "utf8")));
  const shown = (r) => Object.fromEntries((r.samples ?? []).map((s) => [s.sample, { shown: s.shown, views: s.views }]));
  const before = shown(a), after = shown(b);
  const changed = Object.keys(before).filter((k) => JSON.stringify(before[k]) !== JSON.stringify(after[k]))
    .map((k) => ({ sample: k, before: before[k], after: after[k] }));
  console.log(JSON.stringify({ ok: changed.length === 0, rule: "re-running analyze changes nothing a person sees", changed }, null, 2));
  process.exit(changed.length ? 1 : 0);
' "$OUT/3-triaged.json" "$OUT/5-after-analyze.json"

if (( ${#failed[@]} )); then
  echo "verify-truth failed: ${failed[*]}"
  exit 1
fi
echo "verify-truth passed"
