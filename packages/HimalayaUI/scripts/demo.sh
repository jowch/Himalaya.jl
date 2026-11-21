#!/usr/bin/env bash
# Run the real app (Julia backend + Vite dev server) against the synthetic demo
# experiment. Wrapped by `make demo` / `make demo-empty`; see docs/demo.md.
#
#   demo.sh            backend on a persistent demo DB, demo experiment ingested
#   demo.sh --large    same, for the --large tree of experiments in DEMO_DIR
#                      (generated first if missing); every experiment is ingested
#   demo.sh --empty    backend on a FRESH empty DB, nothing ingested, so the
#                      New-experiment funnel can be walked from the start
#
# Env (defaults in brackets):
#   DEMO_DIR       synthetic experiment root      [scratch/demo-experiment]
#   DEMO_DB_DIR    where the demo DB + logs live  [scratch/demo]
#   DEMO_API_PORT  backend port                   [8080]
#   DEMO_UI_PORT   Vite port                      [5173]
#   DEMO_UI_MODE   dev (Vite dev server) or build (vite build, then vite preview) [dev]
# Vite stays in the foreground; Ctrl-C stops both servers.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO"

DEMO_DIR="${DEMO_DIR:-scratch/demo-experiment}"
DEMO_DB_DIR="${DEMO_DB_DIR:-scratch/demo}"
API_PORT="${DEMO_API_PORT:-8080}"
UI_PORT="${DEMO_UI_PORT:-5173}"
JL="julia --project=packages/HimalayaUI"
up() { curl -s -o /dev/null --connect-timeout 1 --max-time 3 "http://127.0.0.1:$API_PORT/api/experiments"; }

empty=0; large=0
case "${1:-}" in
  "") ;;
  --empty) empty=1 ;;
  --large) large=1 ;;
  *) echo "usage: demo.sh [--empty|--large]" >&2; exit 2 ;;
esac

if (( large )); then
  if [[ ! -f "$DEMO_DIR/.demo-large-complete" ]]; then
    echo "No complete large demo at $DEMO_DIR; generating it (make demo-large-data)."
    $JL packages/HimalayaUI/scripts/make_demo_data.jl "$DEMO_DIR" --large
  elif grep -qs ',false,' "$DEMO_DIR"/*/TRUTH.csv; then
    echo "$DEMO_DIR was generated with --no-images (empty TIFFs); regenerate it with make demo-large-data." >&2
    exit 1
  fi
elif [[ ! -f "$DEMO_DIR/TRUTH.csv" ]]; then
  echo "No demo data at $DEMO_DIR; generating it (make demo-data)."
  $JL packages/HimalayaUI/scripts/make_demo_data.jl "$DEMO_DIR"
elif grep -q ',false,' "$DEMO_DIR/TRUTH.csv"; then
  echo "$DEMO_DIR was generated with --no-images (empty TIFFs); run make demo-data to regenerate it." >&2
  exit 1
fi

if (( empty )); then
  DEMO_DB_DIR="$DEMO_DB_DIR-empty"
  rm -rf "$DEMO_DB_DIR"
fi
mkdir -p "$DEMO_DB_DIR"
DB="$(cd "$DEMO_DB_DIR" && pwd)/himalaya.db"

if up; then
  echo "Port $API_PORT is already serving; stop it first (lsof -ti:$API_PORT | xargs -r kill)." >&2
  exit 1
fi

echo "Starting backend on :$API_PORT with DB $DB (log: $DEMO_DB_DIR/backend.log)"
HIMALAYA_DB_PATH="$DB" $JL -e 'using HimalayaUI; main(ARGS)' -- \
  serve --port "$API_PORT" > "$DEMO_DB_DIR/backend.log" 2>&1 &
backend=$!
trap 'kill $backend 2>/dev/null || true' EXIT INT TERM

# Cold start without a sysimage compiles Oxygen + HimalayaUI; allow 5 minutes.
for _ in $(seq 1 300); do
  up && break
  kill -0 $backend 2>/dev/null || { echo "Backend exited; see $DEMO_DB_DIR/backend.log" >&2; exit 1; }
  sleep 1
done
up || { echo "Backend did not answer on :$API_PORT within 5 min; see $DEMO_DB_DIR/backend.log" >&2; exit 1; }

if (( empty )); then
  echo "Empty DB: use '+ New experiment' and point it at $(cd "$DEMO_DIR" && pwd)"
elif (( large )); then
  $JL packages/HimalayaUI/scripts/demo_ingest.jl "$DEMO_DIR"/*/ --port "$API_PORT"
else
  $JL packages/HimalayaUI/scripts/demo_ingest.jl "$DEMO_DIR" --port "$API_PORT"
fi

echo "Open http://127.0.0.1:$UI_PORT/"
cd packages/HimalayaUI/frontend
case "${DEMO_UI_MODE:-dev}" in
  dev) VITE_API_PORT="$API_PORT" npm run dev -- --host 127.0.0.1 --port "$UI_PORT" --strictPort ;;
  build)
    npx vite build
    VITE_API_PORT="$API_PORT" npx vite preview --host 127.0.0.1 --port "$UI_PORT" --strictPort ;;
  *) echo "DEMO_UI_MODE must be dev or build, not ${DEMO_UI_MODE}" >&2; exit 2 ;;
esac
