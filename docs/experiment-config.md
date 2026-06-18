# Experiment Configuration

Every experiment that Himalaya analyses is described by an
`experiment.toml` file living in the experiment directory. This file is
the **source of truth** — it tells Himalaya how to read the lab
notebook, where the data files live, what their naming convention is,
and which beamline parameters were used. The DB is a derived artifact;
the experiment directory is canonical.

## Why a config file

SAXS data lands on disk with conventions that vary by beamline, by
experiment type, and by lab notebook habit. Forcing every scientist to
restructure their lab notebook to match Himalaya's expectations would be
hostile and error-prone. Instead, Himalaya adapts to the lab notebook by
reading a small TOML file that describes its column layout, file
patterns, and beamline parameters. To support a new beamline or
experiment type you write a new TOML — no Julia code changes required.

The config file also captures things the lab notebook doesn't, like the
beamline energy and sample-to-detector flight path, which are needed
for q-calibration and d-spacing calculations downstream.

## Quick start

The canonical workflow:

```bash
# 1. Create a config from the built-in 'simple' template
himalaya config new --type simple --dir /data/ssrl-2025-oct/lipid-a

# 2. Edit experiment.toml to set the experiment name, beamline params,
#    and any column mappings that don't match the simple defaults.
$EDITOR /data/ssrl-2025-oct/lipid-a/experiment.toml

# 3. Initialise — this reads experiment.toml + manifest.csv and
#    creates samples and exposures in the DB.
himalaya init /data/ssrl-2025-oct/lipid-a

# 4. Later, after editing the manifest or fixing the config:
himalaya reingest /data/ssrl-2025-oct/lipid-a
```

The experiment directory is expected to contain at least:

```
/data/ssrl-2025-oct/lipid-a/
  experiment.toml       ← this file
  manifest.csv          ← lab notebook
  data/                 ← raw images
  analysis/
    automatic_analysis/ ← 1D-integrated .dat files
```

## The TOML format

A complete `experiment.toml` has five sections:

```toml
[experiment]
name        = "SSRL-2025-Oct/LipidA"
description = "Lipid A cubic phase screen, October 2025"
manifest    = "manifest.csv"          # relative to this file

[beamline]
energy_kev    = 12.0                  # X-ray beam energy
flight_path_m = 2.5                   # sample-to-detector distance

[manifest]
delimiter      = "\t"                 # "\t" for tab, "," for CSV
skip_rows      = 1                    # rows before the data
header_row     = 0                    # 1-based row of named headers; 0 = positional only

# Each column value is either an integer (1-based position) or a string
# (header name). Named lookup is tried first when header_row > 0.
sample_id      = 1
name           = 2     # stable scientific identifier (e.g. "JC001")
display_name   = 3     # user-facing label (e.g. "DOPC + cholesterol")
filenames      = 9
notes_sample   = 10
notes_exposure = 11

[layout]
data_dir      = "data"
analysis_dir  = "analysis/automatic_analysis"
exposure_type = "simple"

[files]
integration = "{name}.dat"            # how to find a 1D trace from a stem
image       = "{name}.tiff"           # how to find the raw image
```

### `[experiment]`

| Field | Purpose |
|-------|---------|
| `name` | Display name. The convention is `Run/Experiment` (e.g. `"SSRL-2025-Oct/LipidA"`) — Himalaya treats this as an opaque label. If empty, defaults to the directory basename. |
| `description` | Free text, shown in the UI. |
| `manifest` | Path to the lab notebook CSV, relative to `experiment.toml`. May be omitted; Himalaya will register the experiment without samples. |

### `[beamline]`

| Field | Purpose |
|-------|---------|
| `energy_kev` | Beam energy in keV. Used for q-ring radii (q-calibration). |
| `flight_path_m` | Sample-to-detector distance in metres. Used for q-ring radii. |
| `beam_center_x` | Beam center column, in **detector pixels, origin bottom-left** (y-up). |
| `beam_center_y` | Beam center row, in **detector pixels, origin bottom-left** (y-up). |
| `pixel_size_um` | Detector pixel pitch in microns (e.g. Pilatus 172, Eiger 75). |

The Focus page draws physically-correct, beam-centered q-rings only when **all
five** beamline values are present (`energy_kev`, `flight_path_m`,
`beam_center_x`, `beam_center_y`, `pixel_size_um`). Omit any and the rings fall
back to a centered, decorative overlay — no error. Beam center and pixel size
are render-only: they live in the config blob, are surfaced on the experiment
API, and (unlike `energy_kev`/`flight_path_m`) are not mirrored to queryable DB
columns.

`energy_kev` and `flight_path_m` are stored as first-class columns in the DB so
future queries (e.g. "show me all experiments at 12 keV") don't need to parse
the TOML blob.

### `[manifest]`

| Field | Purpose |
|-------|---------|
| `delimiter` | Field separator. `"\t"` (tab) or `","` (comma). |
| `skip_rows` | Number of rows to skip at the top of the file (for preamble like institution headers). |
| `header_row` | 1-based row number containing column headers. `0` means no headers — all column references must be integers. |
| `sample_id`, `name`, `display_name`, `filenames`, `notes_sample`, `notes_exposure` | Column references — see [Column resolution](#column-resolution). |

#### `name` vs `display_name`

The `[manifest].name` column is the **stable scientific identifier** (e.g. `JC001`). It must match
`[A-Za-z0-9._-]+`, be non-empty, and be unique within the experiment. It cannot be edited via the UI;
rename happens through the manifest CSV + reingest.

The `[manifest].display_name` column is the **friendly user-facing label** (e.g. `DOPC + cholesterol`).
It is initialised from the manifest at first ingest and editable via the UI; note that **reingest
refreshes it from the manifest** (alongside `notes`), so a UI edit to `display_name` is overwritten on
the next reingest unless the manifest carries the same value (`cli.jl` `_reingest_inner!`).

Migrating an existing `experiment.toml` from the legacy `label/name` shape to the new
`name/display_name` shape: run `himalaya migrate-toml <experiment-dir>`. Section-aware
regex-anchored substitution; idempotent.

### `[layout]`

| Field | Purpose |
|-------|---------|
| `data_dir` | Where raw images live, relative to the experiment directory. |
| `analysis_dir` | Where 1D-integrated `.dat` files live. |
| `exposure_type` | `"simple"` is the only value currently accepted. Validated at load time against `VALID_EXPOSURE_TYPES` in `config.jl` — unknown values are rejected with a clear error. Reserved for future types like `"raw"`, `"aggregated"`, `"background_subtracted"`. |

### `[files]`

| Field | Purpose |
|-------|---------|
| `integration` | Pattern for the 1D trace file. `{name}` is replaced by the filename stem from the manifest. Resolved relative to `analysis_dir`. |
| `image` | Pattern for the raw image. Resolved relative to `data_dir`. |

Patterns may include subdirectory components — e.g.
`"integrated/{name}.dat"` resolves to `<analysis_dir>/integrated/<stem>.dat`.
Patterns must contain `{name}` and must not be absolute paths or contain
`..` segments — these are rejected at load time with a clear error.

## Column resolution

For each column field in `[manifest]`, the value can be:

- **An integer** — interpreted as a 1-based column index in the row.
- **A string** — interpreted as a header name to look up in the row at
  `header_row`. Falls back to a warning + empty value if the name isn't
  found.

This dual mode lets a lab notebook evolve. A typical positional config:

```toml
[manifest]
header_row     = 0
sample_id      = 1
name           = 2
display_name   = 3
filenames      = 9
```

If the lab adopts named headers later, switch to:

```toml
[manifest]
header_row     = 1
sample_id      = "#"
name           = "ID"
display_name   = "Sample"
filenames      = "Filename(s)"
```

You can mix and match — keep some fields positional and some named —
during a transition.

### Row filtering

Rows whose `sample_id` column doesn't parse as an integer are silently
skipped. This lets scientists add section headers, blank rows, and
free-text annotations to their lab notebook without breaking ingestion.

## Filename-to-exposure association

The manifest's `filenames` column declares which files belong to which
sample. Three forms are supported, all treated as **prefixes against
the filesystem**:

### Range — `JC001-004`

Expanded to `JC001`, `JC002`, `JC003`, `JC004`. For each prefix, the
filesystem is scanned for files matching the integration pattern (e.g.
`JC001*.dat`). All matching files are associated with the sample.

If a file in the range is missing from disk, Himalaya emits a warning
but doesn't fail — scientists sometimes note planned exposures that
didn't complete.

### Single name — `JC001`

The filesystem is scanned for `JC001*.dat`. A single exposure is the
common case, but if the autosampler produced replicates like
`JC001_rep1.dat`, `JC001_rep2.dat`, all of them get associated.

### Bare prefix — `JC`

The filesystem is scanned for `JC*.dat`. All matches are associated
with the sample. This is the "I ran the whole tray and don't want to
list 96 files" case.

In all three cases the **filesystem is the ground truth**. The manifest
declares intent; the disk decides what actually exists.

If the configured `analysis_dir` (or any subdirectory in the integration
pattern) doesn't exist on disk, `resolve_files` emits a `@warn` and
returns no matches — ingestion still succeeds with zero exposures, but
the warning lands in the CLI/server log so the operator can spot a
misconfigured path quickly.

## The read-only contract

Experiment directories are read-only at runtime. Himalaya does not
create, modify, or delete any file inside the experiment directory
during `init`, `analyze`, `reingest`, or `serve`. All write operations
go to the central DB.

The **sole exception** is `himalaya config new --dir <path>`, which
writes `experiment.toml` as a one-time setup step before data
collection. It refuses to overwrite an existing file.

This invariant is enforced by code review and by a regression test:
`test_pipeline.jl` snapshots the directory contents before and after
`cli_init_with_db!` and asserts they're identical.

## Re-ingestion

If the manifest is corrected, the config edited, or new exposures
arrive on disk, run `himalaya reingest <experiment_dir>` to update the
DB. Reingestion is **safe to run repeatedly** and preserves curation
work:

- Existing samples are matched by `(experiment_id, name)` and updated
  in place — both `display_name` and `notes` are **refreshed from the manifest** on reingest.
- Exposures are matched by `(sample_id, filename)`. **Existing exposures
  are never deleted or modified by reingest** — their `accepted` /
  `rejected` status, manual peaks, and analysis results are preserved.
- New exposures discovered on disk are inserted.
- The full TOML blob in `experiments.config` is updated, along with
  `experiment_type`, `energy_kev`, and `flight_path_m`.

The whole operation is wrapped in a SQLite transaction — if anything
goes wrong mid-way, the DB rolls back to its pre-reingest state.

A REST endpoint `POST /api/experiments/:id/reingest` triggers the same
logic from the UI. The response body distinguishes the two reingest
outcomes:

```json
// success — manifest found and ingested
{ "status": "ok", "added_samples": 1, "added_exposures": 12,
  "manifest_path": "/data/.../manifest.csv" }

// success — config updated but no manifest on disk
{ "status": "no_manifest", "added_samples": 0, "added_exposures": 0,
  "manifest_path": "/data/.../manifest.csv" }
```

HTTP status is `200` in both cases — the config-update half of reingest
succeeded, and `status` is the discriminator for callers.

## Built-in templates

Templates live under `packages/HimalayaUI/configs/`. Each `.toml` file
is a starting point that `himalaya config new --type <name>` will copy
into a new experiment directory. The template's content becomes the
scientist's `experiment.toml`, which they then edit.

Currently shipped:

- **`simple.toml`** — encodes the legacy hardcoded behaviour: tab-separated
  Google Sheets export, columns 1/2/3/9/10/11, data in `data/`, .dat
  files in `analysis/automatic_analysis/`, integration files named
  `<stem>.dat`, images named `<stem>.tiff`.

To add a new template, drop a new `.toml` file in `configs/`. It
appears in `himalaya config list` automatically — no code changes
required. Common candidates: `autosampler`, `time_resolved`, beamline-
specific layouts (`bm29`, `id02`, …).

## Backward compatibility

Experiments created before the config system existed have `NULL` in
`experiments.config`. At analysis time, `config_from_db` falls back to
the built-in `simple` template, which exactly matches the old hardcoded
defaults. Existing experiments continue to work without re-initialisation.

If you want to migrate a legacy experiment to a custom config, run
`himalaya config new --dir <exp>` to drop a template, edit it as
needed, and run `himalaya reingest <exp>`. The DB row's `config` blob
gets populated; subsequent `analyze` calls use the new config.

## CLI reference

```
himalaya config list
    List available built-in templates.

himalaya config new [--type TYPE] --dir DIR
    Copy the named template to <DIR>/experiment.toml.
    Defaults to "simple" if --type is omitted.
    Refuses to overwrite an existing experiment.toml.

himalaya init <experiment_dir>
    Read experiment.toml + manifest.csv from the directory and register
    the experiment in the DB. Discovers exposures by filesystem prefix
    scan against the configured integration pattern.

himalaya reingest <experiment_dir>
    Re-read experiment.toml + manifest.csv and update the DB.
    Idempotent on stable input. Preserves curated exposures.

himalaya migrate-toml <experiment_dir>
    Upgrade experiment.toml from the legacy label/name column shape to the new
    name/display_name shape. Section-aware, regex-anchored; idempotent (safe to
    re-run on an already-migrated file). Required once per experiment dir when
    deploying the label/name → name/display_name schema change.
```

## Storage

The full TOML content is stored in `experiments.config` as a TEXT blob,
so the DB is self-contained — `analyze` can resolve file paths even if
the experiment directory is temporarily unmounted. Three additional
columns mirror the most queryable beamline parameters:

| Column | Source |
|--------|--------|
| `experiments.config` | Full `experiment.toml` content |
| `experiments.experiment_type` | `[layout].exposure_type` |
| `experiments.energy_kev` | `[beamline].energy_kev` |
| `experiments.flight_path_m` | `[beamline].flight_path_m` |

The TOML blob is the canonical record; the mirror columns are an
optimisation for queries that don't want to parse TOML.

## Further reading

- [docs/superpowers/specs/2026-04-28-experiment-config-design.md](superpowers/specs/2026-04-28-experiment-config-design.md) — design spec
- [docs/superpowers/plans/2026-04-28-experiment-config.md](superpowers/plans/2026-04-28-experiment-config.md) — implementation plan
- [packages/HimalayaUI/src/config.jl](../packages/HimalayaUI/src/config.jl) — the implementation
- [packages/HimalayaUI/configs/simple.toml](../packages/HimalayaUI/configs/simple.toml) — the default template
