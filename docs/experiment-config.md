# Experiment configuration and ingest

An experiment is a row in the central SQLite DB that points at a beamtime
directory on disk. Its configuration (where the files live, how they are named,
and the beamline geometry) lives on that row. There is no per-experiment
config file and no lab-notebook manifest: an experiment is created through the
New-experiment funnel in the UI, which infers the layout from the directory
tree, lets the user correct it, and then scans the directory.

Experiment directories are read-only data sources. Every write goes to the DB.

## The New-experiment funnel

The funnel is one screen. The user picks a folder, Himalaya lists every setup
under it, and one call per chosen setup creates the experiments. The probes
live in `routes_fs.jl` and never write to the DB; the create and scan live in
`routes_experiments.jl` and `ingest.jl`.

| Step | Route | What it does |
|------|-------|--------------|
| Pick a folder | `GET /api/fs/suggest?prefix=` | Directory autocomplete (at most 20 subdirectories). |
| Preview | `GET /api/fs/preview?path=` | Lists the setups under the folder, each with its layout, counts and geometry (see [Preview](#preview)). `400` when the path is not a directory. |
| Correct a setup | `GET /api/fs/manifest?path=&analysis_dir=&setup_file=&{image,metadata,integration}_pattern=` | Recounts one setup after the user edits its layout in Details (see [File manifest](#file-manifest)). |
| Add | `POST /api/experiments` | One call per chosen setup. Creates the row, returns `202`, and starts the first scan in the background. |

`NewExperimentPage` (`/experiments/new`) shows one row per setup: a checkbox,
an editable name (default `<folder>/<setup>`, e.g. `20261103/setup1`), a
one-line summary, and a Details panel for the folders, patterns and geometry.
A setup that is already an experiment is listed but cannot be chosen.
⌘/Ctrl+Enter adds the chosen setups. One setup opens its workspace straight
away; several go back to the experiments home. The workspace shows the first
scan's progress in its status line (see [Problems](#problems)).

"Manifest" here means the funnel's file manifest: a preview of which images,
metadata files and integration files match the patterns. It is not a lab
notebook.

### Preview

`setup_roots(path)` decides what a setup is. If the picked folder's `data_dir`
(`path/data` if it exists, else `path`) holds `.tif` images, the folder is one
setup. Otherwise each subdirectory whose `data_dir` holds images is a setup, so
picking a beamtime date folder lists `setup1`, `setup2`, and so on, and skips
folders such as `notes`.

For each setup, `resolve_experiment_layout(root)` returns defaults for the user
to correct. It stays cheap because it runs interactively over SMB: a handful of
shallow `readdir` and `isfile` calls, never a deep walk of the raw-image
directory.

| Field | Default |
|-------|---------|
| `name` | `<parent>/<setup>`, the setup folder and the folder above it |
| `data_dir` | `root/data` if it exists, else `root` |
| `analysis_dir` | The directory that holds the integration `.dat` files. Up to 9 image stems, spread evenly through `data_dir`, are each probed against `root/analysis`, its subdirectories, and the subdirectories of `root/analysis/automatic_analysis`. The most common answer wins (ties go to the first seen), so one image whose trace has not been written yet does not decide the layout. Falls back to `root/analysis`, else `nothing`. |
| `image_pattern`, `metadata_pattern`, `integration_pattern` | Set together with `analysis_dir` by the same vote (see [File patterns](#file-patterns)). `layout_detected` is `false` when no stem found an integration file. |
| `setup_file` | The latest `setup_info_*.txt` at the top of `root` or `root/analysis`. |
| `setup_ambiguous` | `true` when no setup file was found, or when setup files exist in both places. |

Each setup also carries the counts a scan would find (`exposures`,
`without_metadata`, `without_trace`), the derived geometry with its sources,
any geometry discrepancies, and `duplicate_of` when an experiment already uses
its `data_dir`.

A new experiment stores the real path of its `data_dir` (symlinks resolved,
trailing slash dropped). The duplicate check resolves only the folder being
picked and compares it against the stored paths as text. It never visits the
folders of existing experiments, so a stored folder on an unreachable network
share cannot stall adding a new one. An experiment created before this rule
under a symlinked path is matched only by that same spelling.

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
spawned task through `run_ingest!` (see [Rescans](#rescans)). It ends with
`ingest_status` set to `complete` or `failed`. Progress goes out on the SSE `curation`
frame with a `kind` of `ingest_started`, `ingest_progress` (with a `stage` of
`discovery` or `analyzing`), and then `ingest_complete` or
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
`experiment_type` are legacy columns that nothing reads (see
[Legacy](#legacy-experimenttoml)). `last_scan_tier` and
`consecutive_empty_ticks` belonged to the deleted rescan timer and are unused.

`PATCH /api/experiments/:id` edits these. The UI no longer calls it: the
Configuration page was removed with the one-screen funnel. A body with no
patchable field is a `400`.

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

1. **SSRL `_tot.dat`.** The beamline writes one `.tif`, `.prp` and `_tot.dat`
   per exposure, including repeated exposures of the same sample at the same
   position. Nothing sums or averages exposures into a `_tot.dat`. If
   `{base}_tot.dat` exists, where `{base}` is the stem with its `_0_NNN`
   suffix removed, the patterns become `{name}_0_NNN.tif`, `{name}_0_NNN.prp`
   and `{name}_tot.dat`, using the suffix of the probed file.
2. **Plain `.dat`.** If `{stem}.dat` exists, the patterns are the defaults.

Every reader of integration files (analysis and the trace route)
resolves the path through `resolve_trace_path` (`pipeline.jl`). It tries the
exact name first, then the name with a trailing `_<rep>_<frame>` removed,
because exposures store the full per-frame stem while the total integration is
named by the acquisition stem.

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

Beam-centred q-rings on the detector image need `energy_kev`,
`flight_path_m`, `beam_center_x`, `beam_center_y` and `pixel_size_um` all set.
No view draws them since the Focus page was deleted (2026-11-20).

## Scan and grouping

`scan_and_group!(db, experiment_id)` is the only code that inserts exposures. It
reads `data_dir`, `analysis_dir` and the patterns from the row, then:

1. **Discovers** exposures with `scan_directory` (`grouping.jl`): one per image
   in `data_dir`, with its metadata file parsed (`parse_prp`, `prp.jl`) and its
   integration file located. If no file matches the image pattern, it returns
   `status = :empty` and writes nothing. When the folder has any metadata files,
   an image whose metadata file has not appeared yet is left for a later pass
   instead of being grouped without a position or time. A folder with no
   metadata files at all is still ingested.
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
   and new files only add rows. A load or sample is created only when it gains
   at least one new exposure, so no pass leaves an empty sample behind. Existing samples and exposures are never
   rewritten, so renames, rejections, peak curation and phase calls survive.
   The one exception is the beam monitors (`transmitted_counts` from
   `I_signal`, `incident_monitor` from `I_0`, `counting_time`).
   `backfill_monitors!` reads them from the stored `prp_path` of every exposure
   still marked `monitors_read = 0`, on every scan request and at the start of
   every pass, so a finished beamtime ingested before the monitors
   were parsed picks them up with no new files. Each file is read once, and a
   stored value is never overwritten.
5. **Analyzes**, outside the write transaction, every exposure that is new,
   never analyzed, whose integration file has appeared since it was last
   tried, or whose stored `analysis_version` is older than `ANALYSIS_VERSION`
   (`pipeline.jl`: the core Himalaya version plus a local counter, bumped when
   the analysis changes). At SSRL BL 4-2 the image and metadata file are written
   when the exposure completes and the integration pipeline writes the
   `_tot.dat` afterwards, so a scan can see an exposure before its trace. Such an
   exposure gets `analysis_error = 'no trace'` and is analyzed by the first pass
   that finds the file. Any other failure is stored in `analysis_error`, with the
   trace file's size and mtime in `analysis_trace_stat`, and the scan moves on.
   It is retried when that file changes (a half-copied trace that is copied
   again), when the version changes, or by a forced pass. A trace rewritten
   after a successful analysis is not picked up.

After a successful pass, `run_ingest!` warms the experiment's thumbnails on a
separate task, after `ingest_complete`, so the workspace opens without waiting
for them. A thumbnail that isn't warm yet renders on first view. A restart
during the warm doesn't resume it: the experiment is already `complete`, so the
next pass with work warms the rest.

The grouping is not edited in the app. A sample split by mistake (for
example a reshoot at a new position) is joined with a bridge in the workspace.

### Problems

`GET /api/experiments/:id/problems` returns `ingest_error`, the message of the
last failed pass (cleared by the next pass that succeeds), and `exposures`:
each exposure still waiting for its trace (`no_trace`) or whose analysis failed
(`failed`), with its sample and reason. The workspace notice summarizes these
(for example "2 exposures waiting for a trace · 1 analysis failed") and its
Details button opens the list; picking a row opens that sample.

## Rescans

Every pass goes through `run_ingest!` (`server.jl`), which holds a per-experiment
lock, so one experiment never has two passes at once; a second caller waits and
then runs its own, usually empty, pass. It sets `ingest_status` to `scanning`,
back-fills beam monitors, runs `scan_and_group!`, and ends `complete` (stamping
`last_scanned_at`) or `failed`. Every pass after the first sends
`phase: "rescan"`.

`ingest_has_work(db, id)` (`ingest.jl`) says whether a pass would change
anything: an image (with its metadata file, when the folder has any) that has no
exposure yet, an exposure whose trace has appeared, one never analyzed, or one
analyzed by an older version, or a failed one whose trace file has changed. It
lists the image, metadata and trace directories once each and stats only
failed exposures, because the scan route runs it inside the request, under the
write lock. `scan_directory` and this check enumerate images through the same
`tif_listing`, so they agree on which exposures exist.

A pass runs:

- **At create**, as above.
- **When the experiment is opened.** The workspace calls
  `POST /api/experiments/:id/scan` once per opened experiment. With no work
  the route answers `200 {status, changed: false}` with the stored
  `ingest_status` and sends no frames; otherwise it answers `202` and runs a
  pass. `{"force": true}` always runs one and also retries every failed
  analysis. The Rescan palette command uses the same route.
- **At server start.** `resume_ingests!` runs a pass for every experiment left
  `scanning` by a restart or crash, and every experiment with work, which
  includes re-analysis after an upgrade that bumps `ANALYSIS_VERSION`.

A pass with an unreachable or empty `data_dir` stops before the analysis step,
so pending traces wait for the data folder as well as the analysis folder.

There is no timer. Data is copied from the beamline and ingested after the
beamtime, so checking on open is enough.

Rescans are additive. A pattern change that stops matching some existing files
does not remove their exposures.

## The read-only contract

Himalaya never creates, modifies or deletes a file inside an experiment
directory while resolving, previewing, scanning, analyzing or serving. The
funnel probes only read, and thumbnails and full images are cached next to
the DB (`<db_dir>/cache/thumb-128/`, `<db_dir>/cache/full-1536/`), not next
to the images. No code path writes into an experiment directory.

This rule is kept by code review. No test checks it for the scan path.

## Legacy: experiment.toml

Before the funnel, each experiment directory carried an `experiment.toml` and a
lab-notebook CSV, and the CLI ingested them. That path, the `upgrade-grouping`
command that retrofitted those experiments to the current grouping, and every
reader of the TOML text kept in `experiments.config` have been deleted. The
column is left in place and never read or sent.

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
```

`-e` takes an id, a name or a path. When it is left out, `analyze` and `show`
use the only registered experiment and fail with a list if there are several.

## Further reading

- [packages/HimalayaUI/src/routes_fs.jl](../packages/HimalayaUI/src/routes_fs.jl): the funnel probes
- [packages/HimalayaUI/src/routes_experiments.jl](../packages/HimalayaUI/src/routes_experiments.jl): create, `PATCH`, scan
- [packages/HimalayaUI/src/ingest.jl](../packages/HimalayaUI/src/ingest.jl): `scan_and_group!`, `ingest_has_work`, `pending_analysis_ids`
- [packages/HimalayaUI/src/grouping.jl](../packages/HimalayaUI/src/grouping.jl): discovery, loads, slots, sample names
- [packages/HimalayaUI/src/geometry.jl](../packages/HimalayaUI/src/geometry.jl): geometry derivation
- [packages/HimalayaUI/src/config.jl](../packages/HimalayaUI/src/config.jl): `integration_pattern_for`
