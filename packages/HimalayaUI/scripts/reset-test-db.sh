#!/usr/bin/env bash
# Reset the live-test backend's DB to a fresh prod-data copy.
#
# post-pr37 tests need a clean DB state per test (each one mutates the
# chosen exposure into a non-reusable state). Called from the spec's
# beforeEach via execFile when RESET_TEST_DB=1.
#
# Env (with defaults that match the local-dev setup at
# /tmp/himalaya-test/):
#   PROD_DB        — source DB to copy (default: $REPO/data/himalaya.db)
#   TEST_DB_PATH   — destination test DB (default: /tmp/himalaya-test/himalaya.db)
#   TEST_BACKEND_PORT — port the backend listens on (default: 8090)
#   TEST_BACKEND_HOST — host the backend binds to (default: 127.0.0.1)
#   HIMALAYAUI_PROJECT — Julia project path (default: $REPO/packages/HimalayaUI)
#   TEST_BACKEND_LOG  — backend log path (default: /tmp/himalaya-test/backend.log)
set -euo pipefail

# Resolve REPO from this script's location (.../packages/HimalayaUI/scripts/).
# When run from a worktree, this resolves to the worktree root, NOT the prod
# /opt/Himalaya.jl checkout. PROD_DB must therefore default to the canonical
# prod location, not $REPO/data/.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"

PROD_DB="${PROD_DB:-/opt/Himalaya.jl/data/himalaya.db}"
TEST_DB_PATH="${TEST_DB_PATH:-/tmp/himalaya-test/himalaya.db}"
TEST_BACKEND_PORT="${TEST_BACKEND_PORT:-8090}"
TEST_BACKEND_HOST="${TEST_BACKEND_HOST:-127.0.0.1}"
HIMALAYAUI_PROJECT="${HIMALAYAUI_PROJECT:-$REPO/packages/HimalayaUI}"
TEST_BACKEND_LOG="${TEST_BACKEND_LOG:-/tmp/himalaya-test/backend.log}"

if [[ ! -f "$PROD_DB" ]]; then
  echo "reset-test-db: PROD_DB not found at $PROD_DB" >&2
  exit 1
fi

# Find any process listening on TEST_BACKEND_PORT and kill it. Tries lsof
# first; falls back to ss. Backend may be slow to release the port — wait.
pid=$(lsof -ti tcp:"$TEST_BACKEND_PORT" 2>/dev/null | head -1 || true)
if [[ -z "${pid}" ]] && command -v ss >/dev/null; then
  pid=$(ss -tlnpH 2>/dev/null | awk -v p=":$TEST_BACKEND_PORT$" '$4 ~ p {match($0,/pid=([0-9]+)/,a); print a[1]; exit}')
fi
if [[ -n "${pid:-}" ]]; then
  kill "$pid" 2>/dev/null || true
  for _ in {1..20}; do
    sleep 0.2
    kill -0 "$pid" 2>/dev/null || break
  done
fi

# Replace the DB.
mkdir -p "$(dirname "$TEST_DB_PATH")"
cp -f "$PROD_DB" "$TEST_DB_PATH"

# Restart backend in the background.
HIMALAYA_DB_PATH="$TEST_DB_PATH" \
  julia --project="$HIMALAYAUI_PROJECT" \
        -e 'using HimalayaUI; main(ARGS)' \
        -- serve --port "$TEST_BACKEND_PORT" --host "$TEST_BACKEND_HOST" \
        > "$TEST_BACKEND_LOG" 2>&1 &
disown

# Wait until the API answers directly.
for _ in {1..60}; do
  if curl -fs "http://$TEST_BACKEND_HOST:$TEST_BACKEND_PORT/api/experiments" \
        > /dev/null 2>&1; then
    # If a Vite proxy is up on TEST_VITE_PORT (default 5180), it likely has
    # stale connections to the killed backend. Hit a couple of endpoints
    # through the proxy to force-reset its connection pool — without this,
    # the FIRST page request after reset takes 30 s+ to respond as Vite
    # retries its dead persistent connections.
    vite_port="${TEST_VITE_PORT:-5180}"
    if curl -fs "http://$TEST_BACKEND_HOST:$vite_port/" > /dev/null 2>&1; then
      curl -fs "http://$TEST_BACKEND_HOST:$vite_port/api/experiments" \
        > /dev/null 2>&1 || true
      curl -fs "http://$TEST_BACKEND_HOST:$vite_port/api/users" \
        > /dev/null 2>&1 || true
    fi
    exit 0
  fi
  sleep 1
done

echo "reset-test-db: backend did not come up within 60 s" >&2
tail -50 "$TEST_BACKEND_LOG" >&2 || true
exit 1
