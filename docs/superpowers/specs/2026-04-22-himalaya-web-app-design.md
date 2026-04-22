# Himalaya Web Application — Design Spec

**Date:** 2026-04-22  
**Status:** Draft — pending user review  
**Scope:** Core web application for semi-automatic SAXS phase indexing

---

## 1. Vision

A single-command web application that ingests a beamtime experiment directory
and manifest, automatically analyzes all SAXS traces, and presents the most
plausible phase assignments for each sample in an interactive interface. The
system should be correct 80% of the time with no user input; the remaining 20%
should be effortless to refine.

**Not in scope for this spec:** monoclinic/tetragonal lattices, cross-experiment
comparison, stacked/waterfall view, summary table page, auto-best-exposure
selection. See `docs/future-feature-ideas.md`.

---

## 2. Architecture

Three layers communicating through well-defined interfaces:

```
Himalaya.jl/                   (git root)
├── Himalaya.jl                (existing core library — minimal changes)
│   peakfinding.jl, index.jl, phase.jl, ...
│
└── packages/
    └── HimalayaUI/            (new sub-package — Julia, Oxygen.jl)
        ├── Project.toml       deps: Himalaya, Oxygen, SQLite, CSV, TOML
        ├── src/
        │   ├── HimalayaUI.jl  entry point + route definitions
        │   ├── db.jl          SQLite schema, migrations, query helpers
        │   ├── manifest.jl    CSV parser → samples + exposures
        │   ├── pipeline.jl    batch analysis orchestration
        │   └── export.jl      CSV / JSON export
        └── frontend/          (TypeScript, Vite, Observable Plot)
            ├── src/
            │   ├── api.ts     typed fetch wrappers for all endpoints
            │   ├── state.ts   client-side state management
            │   └── views/
            │       ├── SampleList.ts
            │       ├── TraceViewer.ts
            │       ├── PropertiesPanel.ts
            │       ├── MillerPlot.ts
            │       └── PhasePanel.ts
            └── dist/          Vite build output, served by Oxygen.jl
```

**Deployment:** single Julia process per server, started with one experiment
path. Oxygen.jl serves the compiled frontend as static files and handles all
`/api/*` routes. Users connect via SSH port-forwarding. No Node.js at runtime.

The `himalaya.db` for a given run always contains exactly one experiment row
(id = 1). The `:id` parameter in experiment API routes is included for
forward-compatibility but is effectively fixed at 1 in v1.

**Development workflow:**
```bash
# Terminal 1 — Julia backend
julia --project=packages/HimalayaUI -e 'using HimalayaUI; serve()'  # Oxygen.jl on :8080

# Terminal 2 — frontend dev server  
cd frontend && npx vite                   # Vite on :5173, proxies /api → :8080
```

**Build for deployment:**
```bash
cd frontend && npx vite build             # outputs to frontend/dist/
julia --project=packages/HimalayaUI -e 'using HimalayaUI; serve()'  # serves dist/ + API from :8080
```

---

## 3. Data Model (SQLite)

One `himalaya.db` per experiment, stored at `<experiment>/himalaya.db`.

```sql
CREATE TABLE users (
    id       INTEGER PRIMARY KEY,
    username TEXT UNIQUE NOT NULL
);

CREATE TABLE experiments (
    id             INTEGER PRIMARY KEY,
    name           TEXT,
    path           TEXT NOT NULL,
    data_dir       TEXT NOT NULL,      -- <experiment>/data/
    analysis_dir   TEXT NOT NULL,      -- <experiment>/analysis/automatic_analysis/
    manifest_path  TEXT,
    created_at     DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE samples (
    id            INTEGER PRIMARY KEY,
    experiment_id INTEGER REFERENCES experiments(id),
    label         TEXT,    -- "D1"
    name          TEXT,    -- "UX1"
    notes         TEXT
    -- 'type' field omitted: Control/Sample distinction handled via tags
);

CREATE TABLE sample_tags (
    id        INTEGER PRIMARY KEY,
    sample_id INTEGER REFERENCES samples(id),
    key       TEXT NOT NULL,    -- "lipid", "peptide", "condition"
    value     TEXT NOT NULL,    -- "DOPC", "melittin", "condensed"
    source    TEXT DEFAULT 'manual'  -- 'manifest' | 'manual'
);

-- A single .dat file, or a derived (averaged / background-subtracted) trace
CREATE TABLE exposures (
    id        INTEGER PRIMARY KEY,
    sample_id INTEGER REFERENCES samples(id),
    filename  TEXT,                        -- NULL for derived exposures
    kind      TEXT DEFAULT 'file',         -- 'file' | 'averaged' | 'background_subtracted'
    selected  BOOLEAN DEFAULT FALSE
);

-- Provenance for derived exposures
CREATE TABLE exposure_sources (
    averaged_exposure_id INTEGER REFERENCES exposures(id),
    source_exposure_id   INTEGER REFERENCES exposures(id),
    role                 TEXT DEFAULT 'signal',  -- 'signal' | 'background'
    PRIMARY KEY (averaged_exposure_id, source_exposure_id)
);

CREATE TABLE exposure_tags (
    id          INTEGER PRIMARY KEY,
    exposure_id INTEGER REFERENCES exposures(id),
    key         TEXT NOT NULL,    -- "flag", "quality"
    value       TEXT NOT NULL,    -- "flare", "poor_signal", "beam_off_center"
    source      TEXT DEFAULT 'manual'
);

CREATE TABLE peaks (
    id          INTEGER PRIMARY KEY,
    exposure_id INTEGER REFERENCES exposures(id),
    q           REAL NOT NULL,
    intensity   REAL,
    prominence  REAL,
    sharpness   REAL,
    source      TEXT DEFAULT 'auto'  -- 'auto' | 'manual'
);

CREATE TABLE indices (
    id          INTEGER PRIMARY KEY,
    exposure_id INTEGER REFERENCES exposures(id),
    phase       TEXT NOT NULL,               -- 'Pn3m', 'Im3m', 'Ia3d', etc.
    basis       REAL NOT NULL,
    score       REAL,
    r_squared   REAL,
    lattice_d   REAL,
    status      TEXT DEFAULT 'candidate'     -- 'candidate' | 'confirmed' | 'excluded' | 'stale'
);

-- Which peaks support which index, and at which ratio position
CREATE TABLE index_peaks (
    index_id       INTEGER REFERENCES indices(id),
    peak_id        INTEGER REFERENCES peaks(id),
    ratio_position INTEGER,   -- index into phaseratios(P) — mirrors SparseVector slot
    residual       REAL,      -- deviation from ideal ratio position
    PRIMARY KEY (index_id, peak_id)
);

-- Append-only audit trail
CREATE TABLE user_actions (
    id          INTEGER PRIMARY KEY,
    user_id     INTEGER REFERENCES users(id),
    timestamp   DATETIME DEFAULT CURRENT_TIMESTAMP,
    action      TEXT,        -- 'add_peak' | 'remove_peak' | 'confirm_index' |
                             --  'exclude_index' | 'select_exposure' | 'add_tag' |
                             --  'remove_tag' | 'analyze'
    entity_type TEXT,        -- 'peak' | 'index' | 'exposure' | 'sample'
    entity_id   INTEGER,
    note        TEXT
);
```

**Key design decisions:**
- `index_peaks.ratio_position` mirrors the `SparseVector` slot in `Index{P}`,
  making Julia↔DB round-trips straightforward.
- `indices.status = 'stale'` is set when a supporting peak is added or removed;
  the frontend prompts re-analysis rather than triggering it automatically.
- `exposure_sources.role` supports both simple averaging ('signal' only) and
  background subtraction ('signal' + 'background') in one table.
- Tags replace boolean flag fields throughout — extensible without schema migrations.

---

## 4. API (Oxygen.jl)

All mutating endpoints record `X-Username` header in `user_actions`.

### Users
```
GET    /api/users                          list all usernames (for dropdown)
POST   /api/users                          register new username
GET    /api/users/:username/actions        audit trail for a user
```

### Experiments
```
POST   /api/experiments                    open/create from {path, manifest_path}
GET    /api/experiments/:id                experiment metadata
PATCH  /api/experiments/:id                update directory paths (beamline config)
POST   /api/experiments/:id/analyze        run full batch analysis
GET    /api/experiments/:id/export         download CSV or JSON summary
```

### Samples
```
GET    /api/experiments/:id/samples        list all samples with tag summary
PATCH  /api/samples/:id                    edit name / notes
POST   /api/samples/:id/tags               add tag {key, value}
DELETE /api/samples/:id/tags/:tag_id       remove tag
```

### Exposures
```
GET    /api/samples/:id/exposures          list exposures for a sample
PATCH  /api/exposures/:id/select           set as selected for analysis
POST   /api/exposures/:id/tags             add tag {key, value}
DELETE /api/exposures/:id/tags/:tag_id     remove tag
POST   /api/exposures/:id/analyze          (re)run analysis on one exposure
```

### Peaks
```
GET    /api/exposures/:id/peaks            all peaks (auto + manual)
POST   /api/exposures/:id/peaks            add manual peak {q}; marks affected indices stale
DELETE /api/peaks/:id                      remove peak; marks affected indices stale
```

### Indices
```
GET    /api/exposures/:id/indices          all candidate indices with status
PATCH  /api/indices/:id                    set status: confirmed | excluded | candidate
```

**No PATCH on peaks** — a peak is fully defined by its q position. Changing a
peak position is modeled as DELETE + POST.

---

## 5. Frontend

### Tech stack
- **TypeScript** compiled by **Vite** (dev + build)
- **Observable Plot** for all charts (log-log traces, Miller-index scatter)
- Vanilla TypeScript for UI components — no framework
- Served from `frontend/dist/` by Oxygen.jl in production

### Page structure
Single page for v1: the per-sample analysis view. Experiment management and
batch operations are handled via the CLI. Additional pages (summary table,
stacked comparison) deferred to future iterations.

### Layout
```
┌──────────────┬───────────────────────────┬─────────────────────┐
│              │      Trace Viewer         │  Miller-Index Plot  │
│              │   (log-log, dominant)     │  (per assignment    │
│ Sample List  │                           │   group; toggleable │
│              │                           │   per candidate)    │
│ [filter/tag] │                           ├─────────────────────┤
│              │                           │  Phase Assignments  │
│ ▶ D1 UX1 ✓  │                           │  ● Pn3m  a=12.3nm  │
│   D2 UX2 ·  │                           │    R²=0.998  ✓ ✗   │
│   D3 UL1 ·  ├───────────────────────────│  ▾ Alternatives (2) │
│   D4 UL2 ·  │     Properties Panel      │    Im3m  a=9.1nm   │
│             │  [Exposures][Peaks]        │    R²=0.71   ✓ ✗   │
│             │  [Tags][Notes]             │                     │
└──────────────┴───────────────────────────┴─────────────────────┘
```

### Sample list (left column)
- One row per sample; status dot: grey (unanalyzed), yellow (candidates only),
  green (at least one confirmed index)
- Filter by tag key or value
- Clicking a row loads that sample into the center and right panels

### Trace viewer (center, top)
- Observable Plot: log-log axes, I(q) with σ error ribbon
- Auto-detected peaks: filled markers; manual peaks: outlined markers,
  distinct color
- Per-confirmed-index overlay: vertical tick marks at predicted q positions,
  colored by phase; hover highlights supporting peaks in right panel
- **Click on trace:** snaps to nearest local maximum within a tolerance window,
  adds a manual peak via `POST /api/exposures/:id/peaks`, marks affected
  indices stale
- **Click on existing peak marker:** removes via `DELETE /api/peaks/:id`,
  marks affected indices stale; stale indices show a "re-analyze" prompt

### Properties panel (center, bottom) — tabbed
- **Exposures:** list of `.dat` files; click to select active; add/remove
  exposure tags; derived exposures show provenance (averaged from / background
  subtracted)
- **Peaks:** table of q, prominence, sharpness, source (auto/manual)
- **Tags:** key-value tag editor for the sample
- **Notes:** free-text notes field, persisted to `samples.notes`

### Miller-index plot (right column, top)
- Observable Plot scatter: x = √(h²+k²), y = observed q
- Linear fit line (slope = lattice parameter d); R² annotation
- Defaults to one combined plot for the confirmed assignment group
- Toggle buttons to show/hide individual candidate indices
- `plot.scale("x").invert(event.offsetX)` used for any future click
  interactions on this plot

### Phase panel (right column, bottom)
- Top-scored assignment(s) shown expanded: phase name, lattice parameter,
  R², confirm (✓) and exclude (✗) buttons
- Alternatives in a collapsible list; same controls
- Confirmed indices highlighted; excluded indices dimmed
- Stale indices show a yellow "re-analyze" badge

### User identification
- Modal on first visit: dropdown of existing usernames from `GET /api/users`,
  plus a "New user" option with a text input
- Selection stored in `localStorage` to pre-select on return visits
- Username sent as `X-Username` header on all mutating requests

---

## 6. CLI

Entry point: `himalaya` (or `julia --project=packages/HimalayaUI -e 'using HimalayaUI; CLI.main()'`)

```bash
himalaya init <experiment_path> --manifest <path>   # create DB, parse manifest
himalaya serve <experiment_path> [--port 8080]       # start web server
himalaya analyze <experiment_path>                   # batch analysis, all samples
himalaya analyze <experiment_path> --sample D1       # re-analyze one sample
himalaya show <experiment_path> --sample D1          # print peaks + top indices
himalaya export <experiment_path> [--format csv|json]  # write summary to stdout
```

`serve` is the bridge between CLI and web UI: it opens the DB, compiles routes,
and starts Oxygen.jl. The web UI takes over from there.

---

## 7. Analysis Pipeline

Triggered by `himalaya analyze` (CLI) or `POST /api/experiments/:id/analyze`
(API). Per-exposure steps:

1. Locate `.dat` file via `analysis_dir` + `exposures.filename`
2. Parse three-column file → `(q, I, σ)` vectors
3. Run `findpeaks(q, I)` → peak positions + prominences + sharpness scores
4. Run `indexpeaks(peaks, proms)` → `Vector{Index}`
5. Auto-group results: partition into non-overlapping assignment groups by
   peak sharing; highest-scored group is the default confirmed assignment
6. Persist peaks → `peaks` table; indices → `indices` + `index_peaks` tables
7. Manual peaks (source = 'manual') are preserved across re-analysis runs;
   auto peaks are replaced

**Re-analysis on manual peak edit:**
- Affected indices (those sharing the edited peak) are marked `status = 'stale'`
- Full re-analysis of that exposure is triggered on user confirmation via
  the "re-analyze" prompt in the UI, or immediately via
  `POST /api/exposures/:id/analyze`

---

## 8. Beamline Configuration

Directory paths are configurable per experiment. Defaults ship for the primary
beamline:

```toml
# config/beamlines/default.toml
[paths]
data_dir     = "data"
analysis_dir = "analysis/automatic_analysis"

[files]
dat_pattern  = "{filename}.dat"    # how to find the .dat file given a filename stem
```

Additional beamline profiles can be added as TOML files. Selected at
`himalaya init` time via `--beamline <name>`.

---

## 9. Manifest Parsing

Input: CSV exported from Google Sheets (or equivalent). Expected columns:

| Column | Field | Notes |
|---|---|---|
| `#` | row index | used to detect section header rows (non-numeric → skip) |
| `Sample` | `samples.label` | "D1" |
| `Name` | `samples.name` | "UX1" |
| `Time(s)` | ignored for now | |
| `Filename(s)` | `exposures.filename` | range "JC001-004" expanded to JC001..JC004 |
| `Notes (Sample)` | `samples.notes` | |
| `Notes (Exposure)` | exposure tag `note` | e.g. "sq", "condensed" |

Section header rows (where `#` is non-numeric or empty) are skipped. Filename
ranges like `JC001-004` and `JC013-JC016` are expanded to individual filenames.
Tags extracted from the manifest are stored with `source = 'manifest'` so
manual additions are not clobbered on re-import.
