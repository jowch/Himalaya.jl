# HimalayaUI

A web application for semi-automatic indexing of SAXS diffraction patterns. Point it at a beamtime experiment directory and a sample manifest, and it will auto-find peaks, index them against known lipid phases, and present the results in a browser for review and refinement.

---

## What you get

- **One command per experiment.** `init` → `analyze` → `serve`.
- **Browser UI** with:
  - Sample list with tag-based filter
  - Log-log trace viewer (I(q) with σ ribbon); click to add/remove peaks
  - Active-group overlay (predicted q positions colored by phase) + hover preview of alternatives
  - Miller-index scatter with linear fit per candidate index
  - Phase panel with Active / Alternatives sections (confirm with `+`, exclude with `−`)
  - Tabbed properties panel: Exposures · Peaks table · Sample tags · Notes
- **Everything persists** in a single `himalaya.db` (SQLite) inside the experiment folder — manifest, auto-picked peaks, manual edits, phase assignments, full audit trail.
- **Multi-user by attribution.** Each browser session identifies as a username; every edit is logged to the `user_actions` audit table.

---

## Requirements

| Tool     | Version |
|----------|---------|
| Julia    | 1.9 or newer |
| Node.js  | 20 or newer (only for the first-time build — deployment runs the pre-built frontend) |
| Browser  | Any modern Chromium/Firefox/Safari |

Runs locally on your workstation or on a lab server over SSH port-forward. No external services.

---

## One-time setup

From the repository root:

```bash
# 1. Resolve Julia dependencies
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.instantiate()'

# 2. Build the frontend (produces packages/HimalayaUI/frontend/dist/)
cd packages/HimalayaUI/frontend
npm install
npm run build
cd -
```

You only need step 2 once per clone (or whenever you pull new frontend changes).

The `himalaya` commands below are shown in their fully-explicit form. If you run them often, create an alias:

```bash
alias himalaya='julia --project=/path/to/Himalaya.jl/packages/HimalayaUI -e "using HimalayaUI; main(ARGS)" --'
# Then: himalaya init ..., himalaya analyze ..., etc.
```

---

## Experiment directory layout

HimalayaUI expects a per-experiment folder laid out roughly like a beamline output:

```
my-experiment/
├── data/                            # raw detector images (not touched by the UI)
├── analysis/
│   └── automatic_analysis/
│       ├── sample-D1-001_tot.dat    # ← the files the UI reads
│       ├── sample-D1-002_tot.dat
│       └── ...
└── himalaya.db                      # created by `himalaya init`
```

Each `_tot.dat` file is a whitespace-separated three-column table of `q  I  σ` (azimuthally integrated trace). No header. The filename stem matches what the manifest references.

The two paths (`data`, `analysis/automatic_analysis`) are the defaults. They can be overridden per-experiment by editing the `experiments` row in the database after `init`, or by a beamline profile (see [Beamline profiles](#beamline-profiles)).

---

## The manifest

The manifest is a **tab-separated** file (a Google Sheets export works directly). The parser expects these columns:

| Column | Field                  | Example                 |
|--------|------------------------|-------------------------|
| 1      | row index              | `1`, `2`, …             |
| 2      | Sample label           | `D1`                    |
| 3      | Sample name            | `UX1`                   |
| 9      | Filename(s)            | `JC001-004`             |
| 10     | Notes (Sample)         | `50% DOPC / 50% DOPE`   |
| 11     | Notes (Exposure)       | `sq`, `condensed`, …    |

Rows whose first column is empty or non-numeric are treated as section headers and skipped, so you can keep human-readable group dividers.

**Filename ranges** like `JC001-004` or `JC013-JC016` are expanded to individual filenames (`JC001`, `JC002`, `JC003`, `JC004`, …). A plain value is used as-is.

---

## CLI reference

All commands accept `--help` for the specific argument list.

### `himalaya init <experiment_path> [--manifest <csv>] [--name <str>]`

Creates `<experiment_path>/himalaya.db`, registers the experiment with its data/analysis paths, and (if `--manifest` is given) imports samples and exposures. Safe to re-run — importing is idempotent against `(sample label, filename)`.

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  init ~/beamtime/2026-04-exp42 \
  --manifest ~/beamtime/2026-04-exp42/manifest.tsv
# → Imported 37 samples from manifest.
# → Initialized experiment #1 at /Users/me/beamtime/2026-04-exp42
```

### `himalaya analyze <experiment_path> [--sample <label>]`

Runs the full pipeline — peak-finding → indexing → auto-grouping → persistence — for every exposure (or only the matching sample). Prints a line per exposure. Idempotent: re-running replaces prior auto-picked peaks and auto groups, but preserves any manual peaks and the user's custom group.

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  analyze ~/beamtime/2026-04-exp42
#   Analyzing D1 / JC001 ... done
#   Analyzing D1 / JC002 ... done
#   Analyzing D2 / JC005 ... done
#   ...

# Or just one sample:
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  analyze ~/beamtime/2026-04-exp42 --sample D1
```

### `himalaya show <experiment_path> --sample <label>`

Prints the stored analysis for one sample — all exposures, their peaks, and their candidate indices. Useful as a quick sanity check without opening the browser.

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  show ~/beamtime/2026-04-exp42 --sample D1

# Exposure: JC001
#   Peaks (5):
#     q=0.0643  prom=1.450  sharp=0.880  [auto]
#     q=0.1114  prom=0.920  sharp=0.640  [auto]
#     ...
#   Indices (3):
#     Pn3m    basis=0.0454  score=1.000  R²=0.9998  d=13.84
#     Im3m    basis=0.0321  score=0.612  R²=0.9823  d= 9.78
#     ...
#   Active group: auto
```

### `himalaya serve <experiment_path> [--port 8080] [--host 127.0.0.1]`

Starts the web server. Blocks until you Ctrl-C. The UI lives at `http://<host>:<port>/` and the JSON API at `/api/*`.

```bash
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  serve ~/beamtime/2026-04-exp42 --port 8080
# → HimalayaUI serving /Users/me/beamtime/2026-04-exp42 on http://127.0.0.1:8080
```

To reach a remote server, forward the port with SSH:

```bash
ssh -L 8080:127.0.0.1:8080 user@lab-workstation
```

Then open `http://localhost:8080/` locally.

---

## Using the web UI

On first visit you'll see a username prompt — enter a new name or pick from existing users. This name is stored in your browser (localStorage) and stamped on every edit you make.

**Left column:** The sample list. Filter by any tag key or value. Click a sample to load it.

**Center top:** The trace viewer (log-log). Auto peaks are filled circles; manual peaks are outlined circles. The phase-colored vertical ticks show the predicted q positions of whichever indices are in the **active group**.
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

## Beamline profiles

Different beamlines lay out their `.dat` output differently. Default paths (`data`, `analysis/automatic_analysis`) are shipped; to override, either:

- Pass different values at init time by editing the row (`UPDATE experiments SET data_dir = ?, analysis_dir = ? WHERE id = 1`), or
- Add a TOML profile under `packages/HimalayaUI/config/beamlines/<name>.toml` and pass `--beamline <name>` at `init` time.

---

## Data model

Everything is stored in `<experiment_path>/himalaya.db` (one SQLite file per experiment). The schema is documented in [the design spec](../../docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md#3-data-model-sqlite). Highlights:

- `peaks` — auto-picked and manually-added peaks. Auto peaks are replaced on re-analysis; manual peaks persist.
- `indices` + `index_peaks` — candidate phase assignments with per-peak ratio positions and residuals. An index's `status` flips to `stale` when the user edits peaks that support it.
- `index_groups` — one `auto` group per exposure always exists; a `custom` group is created on first manual add/remove. Custom wins once present.
- `user_actions` — append-only audit trail keyed by `X-Username`. Use it to answer "who did what, when, and to which exposure".

---

## Troubleshooting

**`no himalaya.db at <path>`** when running `serve`. Run `himalaya init <path> --manifest ...` first.

**`done` prints for every exposure but the UI shows no peaks.** The `_tot.dat` files were probably not found under the registered `analysis_dir`. Check `SELECT data_dir, analysis_dir FROM experiments;` in the SQLite DB and verify the files exist at `<analysis_dir>/<filename>.dat`.

**`SKIP (...)` messages during `analyze`.** The analysis of that exposure raised an error (commonly: missing file, unreadable trace, no peaks). Other exposures continue; the error message indicates the cause.

**Browser stays blank at `/`.** Either the frontend wasn't built (`cd packages/HimalayaUI/frontend && npm run build`) or the server is serving an empty `dist/`. `serve` prints the resolved path on startup.

**Pre-existing analyses "disappear" after re-running `analyze`.** Only **auto** peaks/indices/groups are replaced. Manual peaks and any custom group are preserved. If you see stale-index banners that won't clear, run `analyze` on that exposure (`POST /api/exposures/:id/analyze` from the UI) or globally (`himalaya analyze --sample <label>`).

---

## Developing

For contributing to HimalayaUI itself — test commands, architecture notes, non-obvious gotchas — see [../../CLAUDE.md](../../CLAUDE.md) and the implementation plans under [../../docs/superpowers/plans/](../../docs/superpowers/plans/).
