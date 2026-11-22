# HimalayaUI

A web application for semi-automatic indexing of SAXS diffraction patterns. Point it at a beamtime directory, and HimalayaUI will group the frames into samples, auto-find peaks, index them against known lipid phases, and present the results in a browser for review and refinement.

---

## What you get

- **Create experiments in the browser.** Run `himalaya serve`, choose New experiment, and pick a beamtime directory. HimalayaUI infers the layout, file patterns and geometry for you to confirm, then scans and analyzes the data.
- **A three-pane triage workspace** for reviewing one sample at a time: the sample list with review status, the trace plot (scroll to zoom, click to add or switch off peaks, Log, Linear, Log-log and q²I views), frame thumbnails with drop and restore, and the ranked phase candidates you tick into the sample's call. Comparisons with other samples, buffer subtraction, bridges between samples, labels and a custom index sit in the same screen. ⌘K opens a command palette and `?` lists the keys.
- **Configurable file patterns.** Image, metadata and integration filenames are matched by `{name}` patterns that you confirm when you add a setup, so a different layout needs no code changes.
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
julia --threads=4,1 --project=/path/to/Himalaya.jl/packages/HimalayaUI \
      -e 'using HimalayaUI; main(ARGS)' -- <command> ...
```

---

## Creating an experiment

Start the server (`himalaya serve`), open the UI, and choose **New experiment**. Pick a beamtime date folder (it lists every setup under it, such as `setup1` and `setup2`) or a single setup folder. For each setup HimalayaUI:

1. **Finds** the layout: the raw-image directory (`data/` if present), the directory holding the integration `.dat` files, the file patterns, and the `setup_info_*.txt` calibration file.
2. **Previews** what a scan would find in one line: how many frames, how many still lack a trace, and the flight path from the setup file. Open Details to correct a folder, a pattern or a geometry value.
3. **Scans** when you add it (⌘/Ctrl+Enter): it creates one frame per image, groups frames into loads and samples by time and stage position, and runs peak-finding and indexing on each one. Adding one setup opens its workspace, which shows the scan's progress above the sample list.

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

Runs peak-finding and indexing for every frame that is not dropped (or only the matching sample). Idempotent: re-running replaces prior auto-picked peaks, but preserves manual peaks and phase calls. `--all` runs every registered experiment.

```bash
himalaya analyze -e 1
himalaya analyze -e 1 --sample "HA_85_422 (S01P02)"
```

### `himalaya show [-e <experiment>] --sample <label>`

Prints the stored analysis for one sample: frames, peaks, candidate indices. Useful as a quick sanity check without opening the browser.

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

After a pull, one target does everything:

```bash
cd /opt/Himalaya.jl
git pull
make deploy
```

`make deploy` runs `make setup`, then builds the frontend into `dist.next` and
the sysimage into `build/next/`, and runs `make check-sysimage` on the new
sysimage. Only when all of that passes does it move both into place and run
`sudo systemctl restart himalaya`. It stops at the first step that fails, and
until the swap the running daemon, the CLI and the served UI all stay on the
old version. It takes about 5–10 min. The frontend goes to
`HIMALAYA_FRONTEND_DIST` when `.env` sets it.

It runs `sudo -v` first, so the password prompt comes at the start, not after the
build. sudo's credential cache can expire during the build, so it may ask again
at the restart. To avoid that, allow exactly the restart without a password, e.g.
`<user> ALL=(root) NOPASSWD: /usr/bin/systemctl restart himalaya` in
`/etc/sudoers.d/himalaya`. `HIMALAYA_SERVICE=<unit>` names another systemd unit,
and `HIMALAYA_RESTART=` skips the restart.

Packages and the sysimage go into the shared depot that `bin/himalaya` runs
against: every make target reads `JULIA_DEPOT_PATH` from the same `.env`, so
nothing needs exporting. Schema migrations need nothing extra: `open_db` runs
them on connect, so the restart picks them up.

A frontend-only change doesn't need the full deploy: `make frontend` and a
browser refresh is enough (~30 s).

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

On first visit you'll see a username prompt. Enter a new name or pick an existing one. The name is stored in your browser (localStorage) and stamped on every edit you make.

The experiments home lists every experiment. Opening one goes straight to its workspace on the first unreviewed sample:

- **Left:** the sample list, grouped by setup when there are several, with each sample's call (a phase, "Form factor" or unreviewed) and "N of M reviewed". `/` filters it; `Z` jumps to the next unreviewed sample.
- **Centre:** the trace plot and the frame thumbnails. Scroll zooms q. A click near a peak switches an auto peak off or removes a manual one, and a click elsewhere adds a manual peak. `X` drops the shown frame. `Space` enlarges its detector image. The compare bar above the plot holds up to five comparison samples (`C` adds one), and ⌘/Ctrl+click on a chip subtracts it.
- **Right:** the Call list of ranked candidates (number keys tick them, `0` marks the frame form factor), the indexing plot, and the frame details. `I` opens a custom index.

⌘Z undoes the last edit. The palette (⌘K) holds everything without a key: plot views, bridges, labels and their colours, rescans, settings and theme. The verify skill's [Feature Map](../../.claude/skills/verify-himalaya/references/features/README.md) describes every surface in detail.

---

## Data model

Everything is stored in a single SQLite DB at `default_db_path()` (env-resolved, defaults to `~/.himalaya/himalaya.db`). The schema is defined in [`packages/HimalayaUI/src/db.jl`](src/db.jl). Highlights:

- `experiments` → `setups`: one experiment per beamtime folder, one setup per setup subfolder. A setup holds `data_dir`, `analysis_dir`, the image, metadata and integration patterns, and geometry with a `*_source` per field.
- `loads` → `samples` → `exposures`: the auto-grouping. Loads split by time gaps, samples by stage position within a load. Each exposure row is one frame.
- `auto_peaks` and `peak_curations`: peaks found by analysis, and the user's manual peaks and switched-off auto peaks. Re-analysis replaces auto peaks and keeps curations.
- `indices` + `index_peaks`: candidate phase assignments with per-peak ratio positions and residuals.
- `assignments` + `assignment_members` and `sample_calls`: the phases ticked on a frame, and the sample's call.
- `sample_labels`, `label_colors`, `sample_compares`, `sample_subtractions`, `sample_bridges`: per-sample labels, comparisons, subtraction reference and bridges.
- `user_actions`: append-only event log keyed by `X-Username`. Use it to answer "who did what, when, and to which frame".

---

## Troubleshooting

**Frames or peaks missing after a scan.** If the scan failed, the workspace notice says so; its Details button shows the error and every frame still waiting for a trace or whose analysis failed. To find files that couldn't be paired, call `GET /api/fs/manifest?path=<data_dir>&analysis_dir=<analysis_dir>&image_pattern=…&metadata_pattern=…&integration_pattern=…` with the setup's settings: its `unmatched` list names every image missing its metadata or integration file, with near-miss filenames (wrong case, `.tiff` vs `.tif`). Patterns can't be edited after a setup is added. The only fix is to remove the experiment from the palette and add the folder again with the right patterns. **Removing an experiment deletes every phase call, peak edit, label and bridge you made in it**, though nothing on disk. Do this before triage, not after. Also check that `analysis_dir` points at the directory that actually holds the `.dat` files.

**`SKIP (...)` messages during `analyze`.** The analysis of that frame raised an error (commonly: missing file, unreadable trace, no peaks). Other frames continue; the error message indicates the cause.

**Browser stays blank at `/`.** Either the frontend wasn't built (`cd packages/HimalayaUI/frontend && npm run build`) or `HIMALAYA_FRONTEND_DIST` points somewhere empty. `serve` prints the resolved DB path on startup; check the server logs for `dynamicfiles` mounting.

**Pre-existing analyses "disappear" after re-running `analyze`.** Only auto peaks and candidate indices are replaced. Manual peaks, switched-off peaks, phase calls and custom indices are kept.

---

## Developing

For contributing to HimalayaUI itself — test commands, architecture notes, non-obvious gotchas — see [../../CLAUDE.md](../../CLAUDE.md), [../../docs/experiment-config.md](../../docs/experiment-config.md), and the per-area `AGENTS.md` files (closest to the code you're touching).
