# Frontend dev loop

How to run a HimalayaUI dev session against a live backend. Two scenarios: a single dev session on a clean box (default), or a dev session running side-by-side with live prod / another worktree.

## Default pair (single dev session, no live prod on the box)

Two terminals:

```bash
# Terminal 1 — backend on :8080.
bin/himalaya serve /path/to/experiment --port 8080
# (Or, without sysimage: julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- serve --port 8080)

# Terminal 2 — Vite on :5173, proxying /api/* to :8080.
cd packages/HimalayaUI/frontend && npm run dev -- --host 127.0.0.1
```

Vite reads `VITE_API_PORT` (defaults to `8080`) at config time — see `vite.config.ts`. The `--host 127.0.0.1` flag matters: `playwright.config.ts` expects `http://127.0.0.1:5173`, not `localhost`, and a default-`localhost` Vite will hang Playwright tests for 60 s before failing.

## Side-by-side pair (live prod already on :8080, or two worktrees in parallel)

Use non-default ports + a DB copy. Production data must NEVER be a dev backend's `HIMALAYA_DB_PATH` — copy it first:

```bash
# 1. Snapshot prod DB into a scratch location. Always copy, never share.
mkdir -p /tmp/himalaya-uat
cp /opt/Himalaya.jl/data/himalaya.db /tmp/himalaya-uat/himalaya.db   # adjust source path

# 2. Backend on a free port pointed at the copy.
HIMALAYA_DB_PATH=/tmp/himalaya-uat/himalaya.db \
  julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  serve --port 8091 > /tmp/himalaya-uat/backend.log 2>&1 &

# 3. Vite on a matching free port, with VITE_API_PORT pointing at the dev backend.
#    `VITE_API_PORT=8091` rewrites the proxy target so /api/* → :8091, NOT :8080.
#    Skipping it silently routes dev requests to PROD — easy to lose a half hour
#    debugging "phantom" data, so always set it for the side-by-side case.
cd packages/HimalayaUI/frontend
VITE_API_PORT=8091 npm run dev -- --host 127.0.0.1 --port 5182 \
  > /tmp/himalaya-uat/vite.log 2>&1 &
```

Open `http://127.0.0.1:5182/` in the browser. The browser hits Vite (`:5182`), Vite proxies `/api/*` to `:8091`, the dev backend reads `/tmp/himalaya-uat/himalaya.db`. Prod (`:8080`, prod DB) is untouched.

## Stopping cleanly

Backend and Vite are both detached background processes; kill by listening port. `pkill -f "port <N>"` is exit-144-noisy in some shells but the kill itself succeeds — verify with `ss -tlnp | grep <port>`:

```bash
pkill -f "port 8091" 2>/dev/null   # backend
lsof -ti:5182 | xargs -r kill      # Vite (more reliable than pkill -f for npm wrappers)
ss -tlnp | grep -E "8091|5182" || echo "ports free"
rm -rf /tmp/himalaya-uat            # log + DB copy
```

## Common gotchas

- `localhost` vs `127.0.0.1` — Playwright is strict; bind explicitly.
- A previous dev session left a Vite child on the same port — `lsof -ti:5173 | xargs kill -9` and rerun.
- `VITE_API_PORT` is read at config load time (Vite startup), not per-request — restart Vite if you change the dev backend's port.
- The `--sysimage`-loaded prod backend doesn't pick up new code on file save — kill + restart it after merging a fix you want exercised end-to-end.
