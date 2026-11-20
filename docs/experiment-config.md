# Experiment configuration and ingest

An experiment is a row in the central SQLite DB that points at a beamtime
directory on disk. Its configuration (where the files live, how they are named,
and the beamline geometry) lives on that row. There is no per-experiment
config file and no lab-notebook manifest: an experiment is created through the
New-experiment funnel in the UI, which infers the layout from the directory
tree, lets the user correct it, and then scans the directory.

Experiment directories are read-only data sources. Every write goes to the DB.

## The New-experiment funnel

The funnel is a sequence of read-only filesystem probes followed by one create
call. The probes live in `routes_fs.jl` and never write to the DB; the create
and scan live in `routes_experiments.jl` and `ingest.jl`.

| Step | Route | What it does |
|------|-------|--------------|
| Pick a directory | `GET /api/fs/suggest?prefix=` | Directory autocomplete (at most 20 subdirectories). |
| Gate the pick | `GET /api/fs/validate?path=` | `ok: false` if the path is not a directory or another experiment already uses it as `data_dir`. |
| Resolve | `GET /api/fs/resolve?path=` | Structural auto-discovery of the layout (see [Resolve](#resolve)). |
| Review | `GET /api/fs/manifest?path=&analysis_dir=&setup_file=&{image,metadata,integration}_pattern=` | Previews what a scan would find with the current settings, and derives the geometry (see [File manifest](#file-manifest)). |
| Approve | `POST /api/experiments` | Creates the row, returns `202`, and starts the first scan in the background. |

The frontend holds the in-progress values in a draft store
(`lib/draftExperiment.ts`). `NewExperimentPage` covers the pick, and
`ConfigurationPage` in first-run mode (`/experiments/new/config`) covers resolve,
review and approve. After Approve the user lands on the grouping review page
(`/experiments/:id/grouping`) while the first scan streams progress over SSE.

"Manifest" here means the funnel's file manifest: a preview of which images,
metadata files and integration files match the patterns. It is not a lab
notebook.

### Resolve

`resolve_experiment_layout(root)` returns defaults for the user to correct. It
stays cheap because it runs interactively over SMB: a handful of shallow
`readdir` and `isfile` calls, never a deep walk of the raw-image directory.

| Field | Default |
|-------|---------|
| `name` | `basename(root)` |
| `data_dir` | `root/data` if it exists, else `root` |
| `analysis_dir` | The directory that holds the integration `.dat` files, found by probing the stem of the first `*.tif` in `data_dir` against `root/analysis`, its subdirectories, and the subdirectories of `root/analysis/automatic_analysis`. Falls back to `root/analysis`, else `nothing`. |
| `image_pattern`, `metadata_pattern`, `integration_pattern` | Set together with `analysis_dir` by the same probe (see [File patterns](#file-patterns)). `nothing` when the probe finds no integration file. |
| `setup_file` | The latest `setup_info_*.txt` at the top of `root` or `root/analysis`. |
| `setup_ambiguous` | `true` when no setup file was found, or when setup files exist in both places. |

### File manifest

`GET /api/fs/manifest` answers "what would a scan of these settings find?"
without writing anything. `path` is the `data_dir`. It returns:

- `total`: the number of non-hidden files in `data_dir`.

- `matched`: how many stems match the image, metadata and integration
  patterns. Images and metadata files are matched in `data_dir`, and
  integration files in `analysis_dir` when it is given.
- `unmatched`: every image stem that is missing its metadata or integration
  file, with a `miss` kind (`metadata`, `integration` or `image`) and a `near` filename when an almost-matching file exists (wrong case
  or one extra character in the extension, such as `.PRP` or `.prpx`). Image
  files with a near-miss extension (`.tiff` or `.TIF` against a `.tif` pattern)
  are listed too. These are the usual cause of a scan that finds half the data.
- `geometry` and `discrepancies`: geometry derived from up to 50 metadata files
  plus the `setup_file`, or, when none is given, any `setup_info_*.txt` in
  `data_dir/analysis` (see [Geometry](#geometry)).
- `matched_files`: up to 12 matched image filenames, for display.

### Approve

`POST /api/experiments` takes:

```json
{
  "data_dir": "/data/2026-04/HA",
  "analysis_dir": "/data/2026-04/HA/analysis/automatic_analysis/tot_files",
  "name": "HA",
  "patterns": {
    "image": "{name}_0_001.tif",
    "metadata": "{name}_0_001.prp",
    "integration": "{name}_tot.dat"
  },
  "geometry": {
    "energy_kev": 11.0, "energy_kev_source": "prp",
    "flight_path_m": 1.13, "flight_path_m_source": "setup",
    "beam_center_x": 420.8, "beam_center_x_source": "setup",
    "beam_center_y": 838.8, "beam_center_y_source": "setup",
    "pixel_size_um": 172.0, "pixel_size_um_source": "prp"
  }
}
```

- `data_dir` is required (`path` is still accepted from older callers). It must
  be an existing directory: `400` otherwise, and `409` if another experiment
  already uses it.
- `name` defaults to `basename(data_dir)`. `analysis_dir` defaults to
  `data_dir/analysis` if that exists, else `data_dir`.
- A geometry value with no source is stored with source `user`. Non-numeric
  values and non-string sources are ignored.

The row is created with `ingest_status = 'scanning'` and the route answers
`202 {id, status: "scanning", name, data_dir}`. The first scan then runs on a
spawned task. It ends with `ingest_status` set to `complete` (and the rescan
scheduler armed) or `failed`. Progress goes out on the SSE `curation`
frame with a `kind` of `ingest_started`, `ingest_progress` (with a `stage` of
`discovery`, `analyzing` or `thumbnails`), and then `ingest_complete` or
`ingest_failed`. Frames from a rescan also carry `phase: "rescan"`.

## Configuration fields

All of these are columns on the `experiments` row.

| Column | Set by | Editable later |
|--------|--------|----------------|
| `name` | Create, or `PATCH` | Yes, a plain write |
| `description` | `PATCH` only | Yes, a plain write |
| `data_dir`, `analysis_dir`, `path` | Create | No: `PATCH` returns `400` |
| `image_pattern`, `metadata_pattern`, `integration_pattern` | Create | Yes (see [Rescans](#rescans)) |
| `energy_kev`, `flight_path_m`, `beam_center_x`, `beam_center_y`, `pixel_size_um` | Create, or the scan's fill-only step | Yes: the field's `*_source` becomes `user` |
| `q_units` | `PATCH` only | Yes: its source becomes `user`. Reported as `A-1` with source `default` while `NULL` |

`path` mirrors `data_dir` for new experiments. `manifest_path`, `config` and
`experiment_type` are legacy columns that the funnel leaves `NULL` (see
[Legacy](#legacy-experimenttoml-and-the-config-blob)).

The Configuration page in later-edit mode (`/experiments/:id/config`) edits
these through `PATCH /api/experiments/:id`. A body with no patchable field is a
`400`.

## File patterns

Each pattern contains `{name}`, which stands for an exposure's stem.

| Pattern | Default when the column is `NULL` | Matched in |
|---------|-----------------------------------|------------|
| `image_pattern` | `{name}.tif` | `data_dir` |
| `metadata_pattern` | `{name}.prp` | `data_dir` |
| `integration_pattern` | `{name}.dat` | `analysis_dir` |

The scan enumerates the images, takes each stem, and builds the metadata and
integration paths by putting the stem into those patterns. An exposure exists
for every image, whether or not its metadata or integration file is present.

The resolver recognises two layouts by probing one image stem:

1. **SSRL total integration.** If `{base}_tot.dat` exists, where `{base}` is the
   stem with its `_0_NNN` frame suffix removed, the patterns become
   `{name}_0_NNN.tif`, `{name}_0_NNN.prp` and `{name}_tot.dat`, using the
   frame suffix of the probed file.
2. **Per-frame integration.** If `{stem}.dat` exists, the patterns are the
   defaults.

Every reader of integration files (analysis, the trace route, series traces)
resolves the path through `resolve_trace_path` (`pipeline.jl`). It tries the
exact name first, then the name with a trailing `_<rep>_<frame>` removed,
because older experiments migrated by `upgrade-grouping` store the full
per-frame stem while the total integration is named by the acquisition stem.

## Geometry

Geometry is derived once, during the funnel review, and stored at create with a
source for each field:

| Field | Source | Notes |
|-------|--------|-------|
| `energy_kev` | `prp` | Beam energy from the metadata files. A discrepancy is reported when it varies. |
| `pixel_size_um` | `prp` | Looked up from the detector model in the metadata files. A discrepancy is reported when the model is unknown or there is more than one. |
| `beam_center_x`, `beam_center_y` | `setup` | From the latest `setup_info_*.txt` (AgBe calibration). Detector pixels, origin bottom-left. |
| `flight_path_m` | `setup`, else `prp` | The calibrated mean distance from the setup file, else the nominal pipe length from the metadata files. |
| `q_units` | none | Never derived. Reported as `A-1` with source `default` when unset. |

A field that could not be derived has source `default` and a `NULL` value.

The scan's geometry step is fill-only (`_fill_unset_geometry!`, `ingest.jl`): it
writes a field only while its source is still `default`. It derives from the
scanned metadata files and the setup files it finds by walking up from
`analysis_dir` (`_find_setup_files`), not from the funnel's chosen
`setup_file`. A rescan never
re-derives geometry that create committed or that the user edited.

The Focus page draws beam-centred q-rings only when `energy_kev`,
`flight_path_m`, `beam_center_x`, `beam_center_y` and `pixel_size_um` are all
set. Otherwise it falls back to a centred, decorative overlay.

## Scan and grouping

`scan_and_group!(db, experiment_id)` is the ingest path. (The only other code that
inserts exposures is `regroup_experiment!`, behind `himalaya upgrade-grouping`.) It reads
`data_dir`, `analysis_dir` and the patterns from the row, then:

1. **Discovers** exposures with `scan_directory` (`grouping.jl`): one per image
   in `data_dir`, with its metadata file parsed (`parse_prp`, `prp.jl`) and its
   integration file located. If no file matches the image pattern, it returns
   `status = :empty` and writes nothing.
2. **Fills geometry** that is still unset (see above).
3. **Groups** exposures with `group_into_samples`: it splits them into loads at
   large time gaps, and clusters each load into slots by horizontal position.
   Each slot becomes a sample named `<label> (SNNPMM)`, where the label comes
   from the filename before the `_S<scan>_` token, `S` is the load index and `P`
   is the slot index.
4. **Persists** loads, samples and exposures in one write transaction. Writes
   are insert-only: loads are deduplicated on `(experiment_id, load_index)`,
   samples on `(experiment_id, load_id, slot_index)`, and exposures on
   `(experiment_id, filename)`. A rescan of an unchanged directory is a no-op,
   and new files only add rows. Existing samples and exposures are never
   rewritten, so renames, rejections, peak curation and phase calls survive.
   The one exception is the beam monitors (`transmitted_counts` from
   `I_signal`, `incident_monitor` from `I_0`, `counting_time`): a rescan fills
   any that are still NULL from the metadata file, so exposures ingested before
   the monitors were parsed pick them up. It never overwrites a stored value.
5. **Analyzes** each new exposure with `analyze_exposure!`, outside the write
   transaction, so an analysis failure cannot roll back the ingest. A failure
   is logged and the scan moves on.
6. **Prewarms** thumbnails for the experiment. A failure here never fails the
   scan.

The grouping is a suggestion. The grouping review page and the routes in
`routes_grouping.jl` let the user rename, merge and split samples, and move
exposures between them. Possible mistakes are flagged when the samples are
read: a cross-load merge suggestion when a filename label appears in more than
one load, and a split suggestion when positions inside one sample jump.

## Rescans

There are three ways a scan runs again after create:

- **`POST /api/experiments/:id/scan`** returns `202` and scans in the
  background. With `{"force": true}` it always scans. Otherwise it first runs
  `cheap_change_check`, which compares the number of files in `data_dir` ending
  in the image pattern's suffix with the number of exposure rows, and scans only
  if there are more files. When a scan ran, it puts the scheduler back on
  hourly checks.
- **The rescan scheduler** (`server.jl`) runs the same check on a timer: hourly,
  then daily after 6 quiet hourly checks, then stopped after 3 quiet daily
  checks. A detected change puts it back on hourly, even if that scan then
  fails. The tier and the count of quiet checks are stored on the experiment
  row (`last_scan_tier`, `consecutive_empty_ticks`). The scheduler is armed
  after a successful first scan, and on server start every experiment whose
  ingest is `complete` resumes at its stored tier, so a stopped experiment stays
  stopped until a manual scan finds new files. Time since the last check is not
  stored, so after a restart the first check comes within an hour, then the
  tier's interval applies.
- **A pattern edit.** `PATCH` with a pattern field stores it, and the
  Configuration page then calls the scan route with `force: true` so the new
  pattern applies straight away.

Rescans are additive. A pattern change that stops matching some existing files
does not remove their exposures.

## The read-only contract

Himalaya never creates, modifies or deletes a file inside an experiment
directory while resolving, previewing, scanning, analyzing or serving. The
funnel probes only read, and thumbnails are cached next to the DB
(`<db_dir>/cache/thumb-128/`), not next to the images. No code path writes
into an experiment directory.

This rule is kept by code review. No test checks it for the scan path.

## Legacy: experiment.toml and the config blob

Before the funnel, each experiment directory carried an `experiment.toml` and a
lab-notebook CSV, and the CLI ingested them. That path, including
`himalaya init`, `himalaya reingest`, `POST /api/experiments/:id/reingest` and
`cli_init_with_db!`, has been deleted. What remains:

- **`experiments.config`** holds the TOML text for experiments created that
  way. The funnel leaves it `NULL`.
- **`integration_pattern_for`** (`config.jl`) reads `[files].integration`
  from that blob when the row's `integration_pattern` column is unset, so
  legacy experiments keep finding their traces. Analysis, the trace route and
  series traces resolve the pattern through it.
- **`_beamline_from_config`** (`routes_experiments.jl`) reads beam centre,
  pixel size and `q_units` from the blob for legacy rows with no typed
  geometry columns set. The samples list (`routes_samples.jl`) also falls back
  to the blob for `q_units` whenever the column is `NULL`.
- **`himalaya upgrade-grouping`** retrofits a legacy experiment's samples and
  exposures to the current grouping while keeping its curation. It is a
  dry run unless `--apply` is passed.

Nothing reads an `experiment.toml` file. The commands that wrote one
(`himalaya config new/list`, `himalaya migrate-toml`) and the built-in
`configs/simple.toml` template have been deleted.

## CLI reference

The CLI does not create experiments. It works on experiments already in the
central DB (`HIMALAYA_DB_PATH`, else `~/.himalaya/himalaya.db`).

```
himalaya analyze [-e <experiment> | --all] [-s <sample>]
    Run peak-finding and indexing for every exposure that is not rejected.

himalaya show [-e <experiment>] -s <sample>
    Print the stored peaks, indices and active group for one sample.

himalaya serve [--port N] [--host H]
    Start the web server. Creates an empty DB if none exists, so a new user
    can create the first experiment through the funnel.

himalaya upgrade-grouping [--apply] [-e <experiment>]
    Retrofit legacy experiments to the current auto-grouping. Dry run by
    default.
```

`-e` takes an id, a name or a path. When it is left out, `analyze` and `show`
use the only registered experiment and fail with a list if there are several,
and `upgrade-grouping` runs on every experiment.

## Further reading

- [packages/HimalayaUI/src/routes_fs.jl](../packages/HimalayaUI/src/routes_fs.jl): the funnel probes
- [packages/HimalayaUI/src/routes_experiments.jl](../packages/HimalayaUI/src/routes_experiments.jl): create, `PATCH`, scan
- [packages/HimalayaUI/src/ingest.jl](../packages/HimalayaUI/src/ingest.jl): `scan_and_group!`, `cheap_change_check`, `regroup_experiment!`
- [packages/HimalayaUI/src/grouping.jl](../packages/HimalayaUI/src/grouping.jl): discovery, loads, slots, sample names, flags
- [packages/HimalayaUI/src/geometry.jl](../packages/HimalayaUI/src/geometry.jl): geometry derivation
- [packages/HimalayaUI/src/config.jl](../packages/HimalayaUI/src/config.jl): `integration_pattern_for`
