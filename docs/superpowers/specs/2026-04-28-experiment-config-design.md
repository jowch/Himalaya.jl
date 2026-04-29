# Experiment Config System — Design Spec

**Date:** 2026-04-28  
**Status:** Draft  

## Context

HimalayaUI currently hardcodes all file I/O conventions: data lives in `<exp>/data/`, analysis output goes to `<exp>/analysis/automatic_analysis/`, `.dat` files have no header and three whitespace-separated columns, and the manifest is always a tab-separated Google Sheets export with fields at fixed column positions (1, 2, 3, 9, 10, 11). A `--beamline` CLI flag already exists but is unused.

Different beamlines and experiment types (simple exposures, autosampler, time-resolved) produce data with different folder layouts, filename conventions, and manifest formats. The goal is a config-driven system where an `experiment.toml` file describes how to read a specific experiment's data, replacing all hardcoded assumptions. The config file is the adapter — no code changes are needed to support a new beamline or experiment type.

## Read-Only Experiment Directories

**Experiment directories are read-only at runtime.** Himalaya never creates, modifies, or deletes any file inside an experiment directory during `init`, `analyze`, `reingest`, or `serve`. All write operations go to the central DB (`/opt/himalaya/himalaya.db`).

The **sole exception** is `himalaya config new --dir /path/to/exp`, which writes `experiment.toml` as a one-time setup step before data collection begins. This is explicitly a pre-data operation; once an experiment is active, even this command should not overwrite an existing `experiment.toml`.

This constraint must be verified at every code review touching `pipeline.jl`, `cli.jl`, `config.jl`, or any route handler: no `open(..., "w")`, `write`, `mkdir`, or `rm` call should target a path derived from `experiments.data_dir` or `experiments.analysis_dir`.

## Deployment Model

A single Himalaya instance runs as a daemon on a shared workstation (e.g., in `/opt/himalaya`), with one central SQLite DB. All users access it via SSH port-forwarding. The DB is the primary artifact; experiment data directories are mounted on the workstation but are not required to be present at serve time — only at `init` and `analyze` time.

Experiment data is organized by run (`Run/Experiment` naming convention), but Himalaya treats the experiment as the atomic unit. No `runs` table is needed — the run is encoded in the experiment name string (e.g., `"SSRL-2025-Oct/LipidA"`).

## Experiment Directory Layout

Each experiment directory is self-describing:

```
/data/ssrl-2025-oct/lipid-a/
  experiment.toml       ← config + metadata (source of truth)
  manifest.csv          ← lab notebook / sample sheet
  data/                 ← raw .dat and .tiff files
  analysis/
    automatic_analysis/ ← analysis output
```

`experiment.toml` and `manifest.csv` together are the canonical record. The DB is a derived artifact — a full `reingest` from these files can rebuild it.

## `experiment.toml` Format

Generated from a built-in template via `himalaya config new`. Scientists edit it to match their lab notebook column layout and beamline parameters.

```toml
[experiment]
name        = "SSRL-2025-Oct/LipidA"
description = "Lipid A cubic phase screen, October 2025 SSRL run"
manifest    = "manifest.csv"          # relative to this file's directory

[beamline]
energy_kev    = 12.0                  # X-ray beam energy in keV
flight_path_m = 2.5                   # sample-to-detector distance in meters

[manifest]
delimiter      = "\t"                 # field separator ("\t" or ",")
skip_rows      = 1                    # rows to skip before header/data
header_row     = 2                    # 1-based row of named headers; 0 = no named headers

# Each column value is either an integer (1-based position) or a string (header name).
# Named lookup is attempted first when header_row > 0; falls back to positional on miss.
sample_id      = 1
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
integration = "{name}.dat"            # {name} = filename stem from manifest
image       = "{name}.tiff"           # patterns may include subdirs: "images/{name}.tiff"
```

### Column resolution

For each `[manifest]` field:

1. If the value is a string and `header_row > 0`: look for that name in the header row. If found, use its column index.
2. Otherwise: use the value as a 1-based integer column index.
3. If a named lookup misses: emit a warning and fall back to the built-in `simple` default position.

### File pattern resolution

`{name}` in `[files]` patterns is the filename stem expanded from the manifest. Patterns may include subdirectory components (e.g., `"images/{name}.tiff"`), resolved via `joinpath(data_dir, expanded_pattern)`. Patterns must not be absolute paths or traverse upward (`..`) — validated at `init` time with a clear error.

### Filename-to-exposure association

All filename entries in the manifest are treated as prefixes. The manifest declares which files belong to which sample; the filesystem resolves the actual file list:

- **Range** (`JC001-004`): expand to `["JC001", "JC002", "JC003", "JC004"]`, then for each prefix scan `data_dir` for matching integration files using the `[files].integration` pattern (e.g., if `integration = "{name}.dat"`, scan for `JC001*.dat`).
- **Single entry** (`JC001` or `JC`): scan `data_dir` for all files matching the integration pattern with `{entry}` as prefix.

All matches are sorted by filename and associated with the sample as exposures. A range entry whose files are not found on disk emits a warning but does not fail — scientists sometimes note planned exposures that didn't complete.

## DB Schema Changes

```sql
-- Add to experiments table
ALTER TABLE experiments ADD COLUMN config          TEXT;  -- full experiment.toml content (TOML blob)
ALTER TABLE experiments ADD COLUMN experiment_type TEXT;  -- e.g. "simple"
ALTER TABLE experiments ADD COLUMN energy_kev      REAL;
ALTER TABLE experiments ADD COLUMN flight_path_m   REAL;
```

`config` stores the full `experiment.toml` blob so the DB is self-contained for re-analysis. `energy_kev` and `flight_path_m` are promoted to first-class columns for future q-calibration queries.

## `reingest` Semantics

`himalaya reingest /path/to/exp` re-reads `experiment.toml` and `manifest.csv` and updates the DB:

- Match exposures on `filename` (stable identifier).
- Update sample associations, notes, and config.
- Insert new exposures found on disk.
- **Do not delete** exposures that have `accepted`/`rejected` status or manual peaks — these represent curation work that must not be overwritten by a manifest correction.
- Exposures with no status and no manual peaks that are no longer referenced by the manifest are marked `status = NULL` and left in the DB (not deleted) — they remain visible but unselected.
- Update `experiments.config`, `energy_kev`, `flight_path_m` from the current `experiment.toml`.

## Code Structure

### New file: `packages/HimalayaUI/src/config.jl`

Owns all config logic:

| Function | Purpose |
|----------|---------|
| `load_config(path)` | Read and validate `experiment.toml` from disk; returns `ExperimentConfig` |
| `config_from_db(db, experiment_id)` | Parse stored TOML blob; returns `ExperimentConfig` |
| `resolve_files(config, data_dir, prefix)` | Prefix scan + filesystem lookup; returns sorted filename stems |
| `parse_manifest(config, path)` | Config-driven manifest parser; returns `Vector{ManifestSample}` |

`ExperimentConfig` struct mirrors the TOML structure. `ManifestSample` struct is unchanged — `config.jl` is a pure translation layer.

### New directory: `packages/HimalayaUI/configs/`

Built-in template TOML files shipped with the package:

- `simple.toml` — encodes current hardcoded behavior exactly (backward-compatible default)

Additional templates (e.g., `autosampler.toml`) added as needed, with no code changes required.

### Changes to existing files

| File | Change |
|------|--------|
| `manifest.jl` | Delegates to `config.jl`; `ManifestSample` struct unchanged |
| `pipeline.jl` | `analyze_exposure!` reads path config via `config_from_db` instead of hardcoded strings |
| `cli.jl` | `cli_init` reads `experiment.toml`, calls `load_config`; `--beamline` flag removed; new `config` subcommand added |
| `db.jl` | Schema migration adds `config`, `experiment_type`, `energy_kev`, `flight_path_m` to `experiments` |
| `server.jl` | Expose `reingest` route (`POST /api/experiments/:id/reingest`) |

### New CLI commands

```
himalaya config new [--type simple] --dir /path/to/exp
    Generate experiment.toml from the named built-in template.
    Defaults to "simple" if --type is omitted.

himalaya config list
    List available built-in config types.

himalaya reingest /path/to/exp
    Re-read experiment.toml + manifest.csv for the experiment at this path
    and update the DB. Safe to run multiple times (idempotent on stable data).
```

The existing `himalaya init` command is updated to read `experiment.toml` from the experiment directory (or accept `--config <path>` to point to it explicitly).

## Backward Compatibility

Experiments initialized before this change have no `config` blob in the DB. At `analyze` time, if `config` is NULL, fall back to `simple.toml` defaults (current hardcoded behavior). This means existing experiments continue to work without re-initialization.

## Verification

1. Generate `experiment.toml` with `himalaya config new --dir /tmp/test-exp`
2. Edit column mappings to match a real lab notebook CSV
3. Run `himalaya init /tmp/test-exp` — confirm samples and exposures appear in DB with correct filenames
4. Run `himalaya analyze` on one exposure — confirm `.dat` path is resolved via config, not hardcoded
5. Edit a manifest row, run `himalaya reingest /tmp/test-exp` — confirm sample notes update but curated exposures are preserved
6. After `himalaya init` and `himalaya analyze`, confirm no new files were created inside the experiment directory (only `experiment.toml` written by the prior `config new` step should be present)
7. Run frontend `npm test` and `npm run e2e` — confirm no regressions
8. Run `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'` — all backend tests pass
