# HimalayaUI

A web application for semi-automatic indexing of SAXS diffraction patterns. Drop an `experiment.toml` next to your data and lab notebook, and HimalayaUI will auto-find peaks, index them against known lipid phases, and present the results in a browser for review and refinement.

---

## What you get

- **Five commands per experiment.** `config new` → edit → `init` → `analyze` → `serve`. Plus `reingest` whenever the manifest changes.
- **Browser UI** with:
  - Sample list with tag-based filter
  - Log-log trace viewer (I(q) with σ ribbon); click to add/remove peaks
  - Active-group overlay (predicted q positions colored by phase) + hover preview of alternatives
  - Miller-index scatter with linear fit per candidate index
  - Phase panel with Active / Alternatives sections (confirm with `+`, exclude with `−`)
  - Tabbed properties panel: Exposures · Peaks table · Sample tags · Notes
- **Adapter-driven file I/O.** Different beamlines and experiment types lay out files differently. An `experiment.toml` per experiment describes manifest columns, file patterns, and beamline parameters — no code changes to support a new layout.
- **One central database.** A single SQLite DB at `~/.himalaya/himalaya.db` (override with `HIMALAYA_DB_PATH`) holds every experiment ever registered. Experiment directories are read-only data sources.
- **Re-ingestable.** Edit the manifest or config and run `himalaya reingest` to update the DB. Curation (accepted/rejected exposures, manual peaks) is preserved.
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

## Experiment directory layout

HimalayaUI expects a per-experiment folder roughly like a beamline output, plus an `experiment.toml` and a manifest CSV that you place in the directory:

```
my-experiment/
├── experiment.toml                # config (manifest columns, file patterns, beamline params)
├── manifest.csv                   # lab notebook (TSV / CSV — format described in experiment.toml)
├── data/                          # raw detector images
└── analysis/
    └── automatic_analysis/
        ├── JC001.dat              # ← 1D-integrated traces; what HimalayaUI reads
        ├── JC002.dat
        └── ...
```

Each `.dat` file is a whitespace-separated three-column table of `q  I  σ`. **Experiment directories are read-only at runtime** — HimalayaUI only writes to the central DB. The sole exception is `himalaya config new --dir`, which writes `experiment.toml` once during setup (and refuses to overwrite).

The directory paths and filename patterns above are the **defaults** in `simple.toml`. They're configurable per experiment — see [Configuring an experiment](#configuring-an-experiment).

---

## Configuring an experiment

Every experiment is described by an `experiment.toml`. Generate one from a built-in template:

```bash
himalaya config new --type simple --dir ~/beamtime/2026-04-exp42
# → Created ~/beamtime/2026-04-exp42/experiment.toml from template 'simple'
```

Edit it to set the experiment name, beamline parameters, and column mappings to match your lab notebook:

```toml
[experiment]
name        = "SSRL-2026-Apr/Exp42"
description = "Lipid A cubic phase screen"
manifest    = "manifest.csv"

[beamline]
energy_kev    = 12.0
flight_path_m = 2.5

[manifest]
delimiter      = "\t"          # "\t" for tab, "," for CSV
skip_rows      = 1             # rows before header/data
header_row     = 0             # 1-based row of named headers; 0 = positional only
sample_id      = 1             # column index OR header name (string)
label          = 2
name           = 3
filenames      = 9
notes_sample   = 10
notes_exposure = 11

[layout]
data_dir      = "data"
analysis_dir  = "analysis/automatic_analysis"
exposure_type = "simple"

[files]
integration = "{name}.dat"     # `{name}` = filename stem
image       = "{name}.tiff"    # patterns may include subdirs: "images/{name}.tiff"
```

See [docs/experiment-config.md](../../docs/experiment-config.md) for the full schema, column-resolution semantics, and filename-association rules. To list available templates: `himalaya config list`.

---

## The manifest

The manifest is your lab notebook — typically a tab-separated Google Sheets export. The `[manifest]` section of `experiment.toml` describes how to read it:

- **Positional columns** (integers): the column-1, column-2, etc. position of each field.
- **Named columns** (strings): a header name to look up at the row given by `header_row`. Mix and match — keep some fields positional and some named during a transition.

Rows whose `sample_id` column doesn't parse as an integer are silently skipped, so you can keep human-readable section headers in your notebook.

**Filename ranges** like `JC001-004` or `JC013-JC016` are expanded to individual prefixes (`JC001`, `JC002`, `JC003`, `JC004`). All filename entries — ranges, single names, or bare prefixes — are treated as **prefixes against the filesystem**: each is scanned for matching files using the `[files].integration` pattern. Disk decides what actually exists; missing files emit warnings, not errors.

**Word-boundary matching.** A prefix only matches a filename when the character immediately after the prefix is non-alphanumeric (typically `_`, `.`, or `-`). This prevents the obvious collision where `JC_C04` would otherwise also catch every `JC_C04s_*.dat` file, and the numeric variant where `JC001` would catch `JC0010`. So `JC` alone does **not** match `JC001` — use the range form (`JC001-009`) when you want to enumerate.

---

## CLI reference

All commands accept `--help`. The DB is centralised — every command opens the same database resolved from `HIMALAYA_DB_PATH` or `~/.himalaya/himalaya.db` (created on first use).

### `himalaya config new --type <name> --dir <path>`

Copies the named built-in template to `<path>/experiment.toml`. Refuses to overwrite an existing file. This is the only command that writes to an experiment directory.

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  config new --type simple --dir ~/beamtime/2026-04-exp42
```

### `himalaya config list`

Lists the built-in templates available to `config new`.

### `himalaya init <experiment_path> [--no-analyze]`

Reads `experiment.toml` and the manifest from `<experiment_path>`, registers the experiment/samples/exposures in the central DB, **and runs the full analysis pipeline** (peak-finding + indexing) over every exposure. Discovers exposures by filesystem prefix scan against the integration pattern.

```bash
himalaya init ~/beamtime/2026-04-exp42
# → Imported 37 samples and 148 exposures from manifest.csv.
# → Initialized experiment 'SSRL-2026-Apr/Exp42' (id=1) at /Users/me/beamtime/2026-04-exp42
# → Running analysis (peak-finding + indexing)...
# →   Analyzing D1 / JC_D01_1_S2449 ... done
# →   ...
```

Pass `--no-analyze` to skip the analysis step (run `himalaya analyze -e <id>` later). Useful when you want to inspect what got registered before committing to compute, or when the analysis dir isn't yet populated.

### Identifying an experiment after `init`

Once an experiment is registered, `reingest`, `analyze`, and `show` resolve it
from the central DB rather than from a path argument. They take an
`-e`/`--experiment` flag whose value can be:

- an **id** (`-e 1`)
- a **name** (`-e "SSRL April 2026"`, must match `experiments.name`)
- a **path** (`-e ~/beamtime/2026-04-exp42`, looked up against `experiments.path`)

`-e` is **required** for the write commands (`reingest`, `analyze`) and
**optional** for the read-only `show` (which defaults to the sole registered
experiment).

### `himalaya reingest -e <experiment>`

Re-reads `experiment.toml` + manifest from the experiment's stored path and updates the DB. **Preserves curation** — exposures with `accepted`/`rejected` status or manual peaks are never deleted or modified, only new ones get inserted. Wrapped in a SQLite transaction so partial failures roll back. Safe to run repeatedly.

```bash
himalaya reingest -e 1                           # by id
himalaya reingest -e "SSRL April 2026"           # by name
himalaya reingest -e ~/beamtime/2026-04-exp42    # by path
```

### `himalaya analyze -e <experiment> [--sample <label>]`

Runs the full pipeline — peak-finding → indexing → auto-grouping → persistence — for every exposure (or only the matching sample). Idempotent: re-running replaces prior auto-picked peaks and auto groups, but preserves any manual peaks and the user's custom group.

```bash
himalaya analyze -e 1
himalaya analyze -e 1 --sample D1
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

---

## Environment variables

Deployment is configured through environment variables. See [`.env.example`](.env.example) for the full surface; the load-bearing ones:

| Var | Default | Purpose |
|-----|---------|---------|
| `HIMALAYA_DB_PATH` | `~/.himalaya/himalaya.db` | Central DB location. For workstation deployments, point at `/opt/himalaya/himalaya.db` or similar. |
| `HIMALAYA_CONFIGS_DIR` | bundled `configs/` | Where `config new --type` reads templates from. Drop your own `*.toml` files here for lab-specific defaults. |
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

Everything is stored in a single SQLite DB at `default_db_path()` (env-resolved, defaults to `~/.himalaya/himalaya.db`). Schema overview in [the design spec](../../docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md#3-data-model-sqlite). Highlights:

- `experiments.config` — full `experiment.toml` content as a TEXT blob, so the DB is self-contained for re-analysis even if the experiment directory is unmounted. Mirror columns (`experiment_type`, `energy_kev`, `flight_path_m`) keep beamline params first-class for queries.
- `peaks` — auto-picked and manually-added peaks. Auto peaks are replaced on re-analysis; manual peaks persist.
- `indices` + `index_peaks` — candidate phase assignments with per-peak ratio positions and residuals. An index's `status` flips to `stale` when the user edits peaks that support it.
- `index_groups` — one `auto` group per exposure always exists; a `custom` group is created on first manual add/remove. Custom wins once present.
- `user_actions` — append-only audit trail keyed by `X-Username`. Use it to answer "who did what, when, and to which exposure".

---

## Troubleshooting

**`experiment.toml not found in <path>`** when running `init` or `reingest`. Run `himalaya config new --dir <path>` first, then edit the generated TOML.

**`no database at <path> — run himalaya init first`** when running `serve`. The DB at `default_db_path()` doesn't exist yet. Run `himalaya init <experiment_path>` against your first experiment to create it. Check `HIMALAYA_DB_PATH` if you expected a different location.

**`done` prints for every exposure but the UI shows no peaks.** The `.dat` files were probably not found under the registered `analysis_dir`. Check the experiment row: the `config` column embeds `data_dir` / `analysis_dir` / `[files].integration`. Verify the files actually exist at `<analysis_dir>/<integration_pattern>` and that the manifest filenames are valid prefixes.

**`No integration files found for prefix '<X>' in <dir>`** during `init` or `reingest`. The manifest references files that don't exist on disk under the configured pattern. Either the prefix is wrong, the pattern in `[files].integration` is wrong, or the files truly aren't there yet. The other prefixes still get processed.

**`SKIP (...)` messages during `analyze`.** The analysis of that exposure raised an error (commonly: missing file, unreadable trace, no peaks). Other exposures continue; the error message indicates the cause.

**Browser stays blank at `/`.** Either the frontend wasn't built (`cd packages/HimalayaUI/frontend && npm run build`) or `HIMALAYA_FRONTEND_DIST` points somewhere empty. `serve` prints the resolved DB path on startup; check the server logs for `dynamicfiles` mounting.

**Pre-existing analyses "disappear" after re-running `analyze`.** Only **auto** peaks/indices/groups are replaced. Manual peaks and any custom group are preserved. If you see stale-index banners that won't clear, run `analyze` on that exposure (`POST /api/exposures/:id/analyze` from the UI) or globally (`himalaya analyze --sample <label>`).

**Curation lost after `reingest`.** Shouldn't happen — `reingest` only inserts new exposures, never touches existing ones with status or peaks. If it does, file an issue with the before/after `experiments.config` blob.

---

## Developing

For contributing to HimalayaUI itself — test commands, architecture notes, non-obvious gotchas — see [../../CLAUDE.md](../../CLAUDE.md), [../../docs/experiment-config.md](../../docs/experiment-config.md), and the implementation plans under [../../docs/superpowers/plans/](../../docs/superpowers/plans/).
