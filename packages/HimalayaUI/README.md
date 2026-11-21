# HimalayaUI

A web application for semi-automatic indexing of SAXS diffraction patterns. Point it at a beamtime directory, and HimalayaUI will group the exposures into samples, auto-find peaks, index them against known lipid phases, and present the results in a browser for review and refinement.

---

## What you get

- **Create experiments in the browser.** Run `himalaya serve`, choose New experiment, and pick a beamtime directory. HimalayaUI infers the layout, file patterns and geometry for you to confirm, then scans and analyzes the data.
- **Browser UI** with:
  - Sample list with tag-based filter
  - Log-log trace viewer (I(q) with σ ribbon); click to add/remove peaks
  - Active-group overlay (predicted q positions colored by phase) + hover preview of alternatives
  - Miller-index scatter with linear fit per candidate index
  - Phase panel with Active / Alternatives sections (confirm with `+`, exclude with `−`)
  - Tabbed properties panel: Exposures · Peaks table · Sample tags · Notes
- **Configurable file patterns.** Image, metadata and integration filenames are matched by `{name}` patterns that you can correct per experiment, so a different layout needs no code changes.
- **One central database.** A single SQLite DB at `~/.himalaya/himalaya.db` (override with `HIMALAYA_DB_PATH`) holds every experiment ever registered. Experiment directories are read-only data sources.
- **Additive rescans.** New files, and traces that land after their images, are picked up when the experiment is opened, on Rescan, and at server start. Rescans only add rows, so curation (rejections, renames, manual peaks, phase calls) is preserved.
- **Multi-user by attribution.** Each browser session identifies as a username; every edit is logged to the `user_actions` audit table.

---

## Requirements

| Tool     | Version |
|----------|---------|
| Julia    | 1.9 or newer |
| Node.js  | 20 or newer (only for the first-time build — deployment runs the pre-built frontend) |
| Browser  | Any modern Chromium/Firefox/Safari |

Runs locally on your workstation or on a lab server over SSH port-forward. No external services.

`bin/himalaya` and the `Makefile` require bash and GNU Make — Linux and macOS only.

---

## One-time setup

From the repository root:

```bash
# 1. Register the lab's private package registry (one-time per Julia depot —
#    Himalaya.jl resolves through it). Skip if it's already registered.
julia -e 'using Pkg; Pkg.Registry.add(Pkg.RegistrySpec(url="https://github.com/Wong-Lab/Registry.jl"))'

# 2. Resolve Julia dependencies
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.instantiate()'

# 3. Build the frontend (produces packages/HimalayaUI/frontend/dist/)
make frontend
```

You only need step 3 once per clone (or whenever you pull new frontend changes).
Step 1 is per-depot, not per-clone — if you set `JULIA_DEPOT_PATH` (e.g. for a
shared multi-user deploy, see [.env.example](.env.example)), run step 1 with
that env var set so the registry lands in the shared depot.

**Optional: build a sysimage for fast startup (~15× speedup, ~5 min one-time cost)**

```bash
make sysimage          # → build/himalaya.so  (resolves build deps automatically)
make check-sysimage    # verify it matches current Julia version
```

Rebuild the sysimage after upgrading Julia — `make check-sysimage` will tell you when it's stale.

**`himalaya` command**

The repo ships a wrapper at `bin/himalaya` that automatically uses the sysimage when present and falls back to plain Julia otherwise:

```bash
# Add to PATH (symlink onto a directory already in PATH)
sudo ln -s /path/to/Himalaya.jl/bin/himalaya /usr/local/bin/himalaya

# Or add the repo's bin/ to your PATH
export PATH="/path/to/Himalaya.jl/bin:$PATH"
```

The examples below use `himalaya` assuming it is on your PATH. Substitute the explicit form if preferred:

```bash
julia --project=/path/to/Himalaya.jl/packages/HimalayaUI \
      -e 'using HimalayaUI; main(ARGS)' -- <command> ...
```

---

## Creating an experiment

Start the server (`himalaya serve`), open the UI, and choose **New experiment**. Pick a beamtime date folder (it lists every setup under it, such as `setup1` and `setup2`) or a single setup folder. For each setup HimalayaUI:

1. **Finds** the layout: the raw-image directory (`data/` if present), the directory holding the integration `.dat` files, the file patterns, and the `setup_info_*.txt` calibration file.
2. **Previews** what a scan would find in one line: how many exposures, how many still lack a trace, and the flight path from the setup file. Open Details to correct a folder, a pattern or a geometry value.
3. **Scans** when you add it (⌘/Ctrl+Enter): it creates one exposure per image, groups exposures into loads and samples by time and stage position, and runs peak-finding and indexing on each one. Adding one setup opens its workspace, which shows the scan's progress above the sample list.

A typical SSRL directory looks like this:

```
my-experiment/
├── data/                                  # raw detector images (.tif) + metadata (.prp)
└── analysis/
    ├── setup_info_20260412_101500.txt     # AgBe calibration: beam centre, distance
    └── automatic_analysis/
        └── tot_files/
            └── HA_85_422_S2404_tot.dat    # 1D-integrated trace: q  I  σ
```

**Experiment directories are read-only at runtime.** HimalayaUI writes only to the central DB, and caches thumbnails next to it.

See [docs/experiment-config.md](../../docs/experiment-config.md) for the funnel, the configuration fields, file patterns, geometry, grouping and rescans.

---

## CLI reference

All commands accept `--help`. The DB is centralised — every command opens the same database resolved from `HIMALAYA_DB_PATH` or `~/.himalaya/himalaya.db` (created on first use).

The CLI does not create experiments; use the New-experiment funnel in the UI.

### Identifying an experiment

`analyze` and `show` take an `-e`/`--experiment` flag whose value can be:

- an **id** (`-e 1`)
- a **name** (`-e "SSRL April 2026"`, must match `experiments.name`)
- a **path** (`-e ~/beamtime/2026-04-exp42`, looked up against `experiments.path`)

Without `-e`, `analyze` and `show` use the only registered experiment and fail with a list if there are several.

### `himalaya analyze [-e <experiment> | --all] [--sample <name>]`

Runs peak-finding and indexing for every exposure that is not rejected (or only the matching sample). Idempotent: re-running replaces prior auto-picked peaks, but preserves manual peaks and phase calls. `--all` runs every registered experiment.

```bash
himalaya analyze -e 1
himalaya analyze -e 1 --sample "HA_85_422 (S01P02)"
```

### `himalaya show [-e <experiment>] --sample <label>`

Prints the stored analysis for one sample — exposures, peaks, candidate indices. Useful as a quick sanity check without opening the browser.

### `himalaya serve [--port 8080] [--host 127.0.0.1]`

Starts the web server. Serves all experiments in the central DB; no per-experiment specifier is needed. Blocks until you Ctrl-C. The UI lives at `http://<host>:<port>/` and the JSON API at `/api/*`.

```bash
himalaya serve --port 8080
# → HimalayaUI serving DB at /Users/me/.himalaya/himalaya.db on http://127.0.0.1:8080
```

To reach a remote server, forward the port with SSH:

```bash
ssh -L 8080:127.0.0.1:8080 user@lab-workstation
```

Then open `http://localhost:8080/` locally.

If no DB exists at the resolved path, `serve` creates an empty one and warns, so a new user can create the first experiment in the UI.


---

## Environment variables

Deployment is configured through environment variables. See [`.env.example`](.env.example) for the full surface; the load-bearing ones:

| Var | Default | Purpose |
|-----|---------|---------|
| `HIMALAYA_DB_PATH` | `~/.himalaya/himalaya.db` | Central DB location. For workstation deployments, point at `/opt/himalaya/himalaya.db` or similar. |
| `HIMALAYA_HOST` | `127.0.0.1` | `serve` bind address (override only if behind an auth-aware reverse proxy). |
| `HIMALAYA_PORT` | `8080` | `serve` bind port. |
| `HIMALAYA_FRONTEND_DIST` | bundled `frontend/dist/` | Path to the built frontend (for `/opt`-style deployments where build artefacts ship separately). |

Julia doesn't auto-load `.env` files. Use `direnv`, source them in your shell, or pass them inline:

```bash
HIMALAYA_DB_PATH=/opt/himalaya/himalaya.db himalaya serve
```

---

## Updating a deployment

The work needed depends on what changed:

| Change type | Steps | Cost |
|---|---|---|
| Frontend only (TS/CSS/React) | `git pull` → `make frontend` → browser refresh | ~30 s |
| Backend Julia, no new deps | `git pull` → `make sysimage` → `sudo systemctl restart himalaya` | ~5–10 min |
| Backend Julia, no sysimage | `git pull` → `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.precompile()'` → `sudo systemctl restart himalaya` | ~2 min, before the restart |
| Backend with new Julia deps | as above + `Pkg.instantiate()` against the shared depot **before** `make sysimage` | ~5–10 min |
| Schema migration | none extra — `open_db` runs migrations on connect, daemon restart picks them up | — |
| Julia version bump (juliaup) | rebuild sysimage (`make check-sysimage` flags the mismatch) | ~5–10 min |

For multi-user deploys, set `JULIA_DEPOT_PATH` for the build steps so packages
and the sysimage write into the shared depot:

```bash
cd /opt/Himalaya.jl
git pull
export JULIA_DEPOT_PATH=/opt/Himalaya.jl/.julia        # or whatever .env has
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.instantiate()'   # only if Project.toml changed
make frontend                                          # only if frontend/ changed
make sysimage                                          # only if any Julia code changed
sudo systemctl restart himalaya                        # only if you rebuilt the sysimage
```

A few non-obvious things to know:

- **Frontend updates are hot.** Oxygen's static-file mount reads `dist/` per
  request, so curators see new UI on browser refresh — no daemon restart.
- **Sysimage swap doesn't kill the running daemon.** systemd keeps the
  already-loaded sysimage in RAM; you can `make sysimage` (which overwrites
  `build/himalaya.so` in place) without downtime. CLI invocations starting
  during the rebuild will see the new file once it's fully written.
- **Precompile before you restart.** Without a sysimage, the first start after a
  backend change precompiles HimalayaUI, including a warm-up that drives every
  workspace route once (~100 s). Run `Pkg.precompile()` while the old daemon is
  still serving, and the restart then answers in ~8 s. With a sysimage, `make
  sysimage` does this, and `bin/himalaya` keeps using `build/himalaya.so` until
  you rebuild it, so a pull without `make sysimage` keeps serving the old code.
  Don't copy a dev `LocalPreferences.toml` (written by `make dev-setup`) to the
  server: it turns the warm-up off.
- **CLI users do nothing.** The wrapper resolves the sysimage on each
  invocation; whoever runs `himalaya …` next gets the new code.
- **DB migrations are idempotent.** `create_schema!` is `CREATE TABLE IF NOT
  EXISTS`; column-add and AUTOINCREMENT migrations only fire when needed.

---

## Using the web UI

On first visit you'll see a username prompt — enter a new name or pick from existing users. This name is stored in your browser (localStorage) and stamped on every edit you make.

**Left column:** The sample list. Filter by any tag key or value. Click a sample to load it.

**Center top:** The trace viewer (log-log). Auto peaks are filled triangles in ice-blue; manual peaks are magenta. A **track row** above the plot carries phase-colored ticks at every predicted-q position for the **active group** — this is the legend by position. Inside the plot area, predicted-q vertical lines stay neutral gray so they don't compete with the trace data; **hovering an index card** lights up that index's ticks in phase color (in both the track and the plot) and dims the others. A **legend strip** under the plot labels the peak triangle types and adds a "predicted {Phase}" entry while you hover.
- **Click on the trace** → adds a manual peak, snapped to the nearest local maximum.
- **Click on a peak marker** → removes that peak.
- Any peak edit marks the affected indices stale; a banner above the trace offers a "Re-analyze" button.

**Center bottom:** Tabbed properties panel:
- **Exposures** — one row per `.dat` file; click to load it into the viewer.
- **Peaks** — table with q, prominence, sharpness, source.
- **Tags** — add or remove sample-level tags (key/value pairs).
- **Notes** — free-text notes; saved on blur.

**Right top:** Miller-index scatter. x = √(h²+k²+l²) (normalized), y = observed q. Each index contributes one point per matched peak; the linear fit's slope is the lattice parameter d.

**Right bottom:** Phase panel.
- **Active group** — the leading assignment. The first time you click `+` or `−`, a "custom" group is cloned from the auto group and becomes active; the auto group is preserved as an alternative.
- **Alternatives** — other candidate indices ranked by score. Hovering an alternative speculatively draws its predicted positions on the trace. Click `+` to add to the active group.

---

## Data model

Everything is stored in a single SQLite DB at `default_db_path()` (env-resolved, defaults to `~/.himalaya/himalaya.db`). The schema is defined in [`packages/HimalayaUI/src/db.jl`](src/db.jl). Highlights:

- `experiments` — one row per experiment: `data_dir`, `analysis_dir`, the image/metadata/integration patterns, and geometry with a `*_source` per field. `config` holds the TOML of legacy experiments and is `NULL` for new ones.
- `loads` → `samples` → `exposures` — the auto-grouping: loads split by time gaps, samples by stage position within a load.
- `peaks` — auto-picked and manually-added peaks. Auto peaks are replaced on re-analysis; manual peaks persist.
- `indices` + `index_peaks` — candidate phase assignments with per-peak ratio positions and residuals. An index's `status` flips to `stale` when the user edits peaks that support it.
- `index_groups` — one `auto` group per exposure always exists; a `custom` group is created on first manual add/remove. Custom wins once present.
- `user_actions` — append-only audit trail keyed by `X-Username`. Use it to answer "who did what, when, and to which exposure".

---

## Troubleshooting

**Exposures or peaks missing after a scan.** If the scan failed, the workspace notice says so; its Details button shows the error and every exposure still waiting for a trace or whose analysis failed. To find files that couldn't be paired, call `GET /api/fs/manifest?path=<data_dir>&analysis_dir=<analysis_dir>&image_pattern=…&metadata_pattern=…&integration_pattern=…` with the experiment's settings: its `unmatched` list names every image missing its metadata or integration file, with near-miss filenames (wrong case, `.tiff` vs `.tif`). Fix the pattern with `PATCH /api/experiments/:id`; a pattern edit forces a rescan. Also check that `analysis_dir` points at the directory that actually holds the `.dat` files.

**`SKIP (...)` messages during `analyze`.** The analysis of that exposure raised an error (commonly: missing file, unreadable trace, no peaks). Other exposures continue; the error message indicates the cause.

**Browser stays blank at `/`.** Either the frontend wasn't built (`cd packages/HimalayaUI/frontend && npm run build`) or `HIMALAYA_FRONTEND_DIST` points somewhere empty. `serve` prints the resolved DB path on startup; check the server logs for `dynamicfiles` mounting.

**Pre-existing analyses "disappear" after re-running `analyze`.** Only **auto** peaks/indices/groups are replaced. Manual peaks and any custom group are preserved. If you see stale-index banners that won't clear, run `analyze` on that exposure (`POST /api/exposures/:id/analyze` from the UI) or globally (`himalaya analyze --sample <label>`).

---

## Developing

For contributing to HimalayaUI itself — test commands, architecture notes, non-obvious gotchas — see [../../CLAUDE.md](../../CLAUDE.md), [../../docs/experiment-config.md](../../docs/experiment-config.md), and the per-area `AGENTS.md` files (closest to the code you're touching).
