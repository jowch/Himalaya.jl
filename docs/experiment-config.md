# Experiment configuration and ingest

An experiment is a row in the central SQLite DB that points at a beamtime
directory on disk. Each setup subfolder of that directory is a row in the
`setups` table under the experiment, and the setup holds the configuration:
where the files live, how they are named, and the beamline geometry. Loads,
samples and exposures carry the `setup_id` they were scanned from. There is no
per-experiment config file and no lab-notebook manifest: an experiment is
created through the New-experiment funnel in the UI, which infers the layout
from the directory tree, lets the user correct it, and then scans the directory.

Experiment directories are read-only data sources. Every write goes to the DB.

## The New-experiment funnel

The funnel is one screen. The user picks a folder, Himalaya lists every setup
under it, and one call creates one experiment for the folder with all of its
setups. The probes
live in `routes_fs.jl` and never write to the DB; the create and scan live in
`routes_experiments.jl` and `ingest.jl`.

| Step | Route | What it does |
|------|-------|--------------|
| Browse | `GET /api/fs/list?path=` | One folder's subfolders, natural order, hidden ones left out, `~` expanded, plus its `parent` (`null` at `/`) and `inside_setup`, true when the folder is a setup or holds images itself. With no path it opens beside the most recently added beamtime, or in the home folder. One directory read, typed from the directory entries, so a big or slow folder never costs a stat per entry. `400` for a relative path or a missing folder. |
| Label a folder | `GET /api/fs/describe?path=` | Marks one folder `beamtime` (with its setup count) or `setup`, and names the experiment it was added as. The browser asks for each listed folder in turn, after the list is shown, and asks for no more once you leave the folder (the one already asked still finishes on the server). Inside a setup it asks for none, so `data` and `analysis` are never read. Finding a setup stops at the first `.tif` in its data folder, but a setup whose images sit directly in it (no `data/`) is read in full once, to look for setup subfolders. |
| Preview | `GET /api/fs/preview?path=` | Lists the setups under the folder, each with its layout, counts and geometry (see [Preview](#preview)). `400` when the path is not a directory. |
| Correct a setup | `GET /api/fs/manifest?path=&analysis_dir=&setup_file=&{image,metadata,integration}_pattern=` | Recounts one setup after the user edits its layout in Details (see [File manifest](#file-manifest)). |
| Add | `POST /api/experiments` | One call for the folder. Creates the experiment and its setups, returns `202`, and starts the first scan in the background. |

The folder browser (`components/FolderPicker.tsx`, drawn by `ui/FolderBrowser.tsx`)
follows the list keys: ↑↓ or Tab move, Enter adds a beamtime or setup (or opens a
plain folder), → opens, ← or Esc go up, typing filters by name and `/` opens a
field for pasting a path. The Data and Integration folder fields under Details
open the same browser, where Enter uses the highlighted folder and Esc closes it.
A start folder that no longer exists opens the nearest folder above it.

Every folder read in these routes and in the scan goes through `dir_entries.jl`
(`list_names`, `folder_names`, `any_name`). Julia's own `readdir` doesn't let
garbage collection run while the directory is read, so on a slow network share
one listing can stall every thread that allocates, ingest and every request
included. That is the likely cause of the slow folder suggestions during an
ingest: measured on local disk, not yet on the share. The `dir_entries.jl`
readers make the libuv call `gc_safe`.

`NewExperimentPage` (`/experiments/new`) shows a name field (default: the
picked folder's name) and one row per setup: a checkbox, the setup's folder
name, a one-line summary, and a Details panel for the folders, patterns and
geometry. A setup that is already in an experiment is listed but cannot be
chosen, and is not sent. Unticked setups are sent with `included: false`, so the
experiment remembers them as hidden and never scans them. ⌘/Ctrl+Enter adds the
experiment and opens its workspace. The workspace shows the first
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
The experiments home shows that real path, not the spelling that was picked,
and rescans read it. If a share is reached through a link whose target changes
between mounts, remove the experiment and add the folder again.

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

`POST /api/experiments` takes the picked folder and every setup to store:

```json
{
  "name": "20260426",
  "path": "/data/20260426",
  "setups": [
    {
      "root": "/data/20260426/setup1",
      "data_dir": "/data/20260426/setup1/data",
      "analysis_dir": "/data/20260426/setup1/analysis/automatic_analysis/tot_files",
      "included": true,
      "patterns": {
        "image": "{name}_0_001.tif",
        "metadata": "{name}_0_001.prp",
        "integration": "{name}_tot.dat"
      },
      "geometry": {
        "energy_kev": 11.0, "energy_kev_source": "prp",
        "flight_path_m": 1.13, "flight_path_m_source": "setup"
      }
    },
    { "root": "/data/20260426/setup2", "data_dir": "/data/20260426/setup2/data", "included": false }
  ]
}
```

- `path` must be an existing directory and is stored as its real path. `name`
  defaults to the folder's name; surrounding spaces are trimmed.
- `setups` must be a non-empty list with at least one `included` setup
  (`included` defaults to `true`). Setups are stored in the order sent, as
  `position` 0, 1, … Each `data_dir` must be an existing directory: `400`
  otherwise, `400` when two setups share one, and `409` if another experiment
  already uses it. `root` defaults to the folder that holds `data_dir`; it must
  exist, hold `data_dir` and sit inside `path` (`400` otherwise). The setup's
  `folder` is the root's name. A setup that is not an object, has no string
  `data_dir`, or has a non-string `root`/`analysis_dir`, a non-boolean
  `included`, or `patterns`/`geometry` that are not objects is a `400`, and so
  is a non-string `name`.
- `path` may not already be another experiment's `path`: that is a `409` with
  `experiment_id` and `name`, and the funnel then offers to open that
  experiment instead. `GET /api/fs/preview` reports the same as `existing`, so
  the funnel says so before anything is sent. New setup folders in a beamtime
  that is already added come in through a rescan.
- `analysis_dir` defaults to `data_dir/analysis` if that exists, else `data_dir`.
- A geometry value with no source is stored with source `user`. Non-numeric
  values and non-string sources are ignored.

A body without `setups` is the older one-setup form: `data_dir` (or `path`),
`name` (default `basename(data_dir)`), `analysis_dir`, `patterns` and
`geometry` at the top level make an experiment with one setup.

The row is created with `ingest_status = 'scanning'` and the route answers
`202 {id, status: "scanning", name, data_dir}`. The first scan then runs on a
spawned task through `run_ingest!` (see [Rescans](#rescans)). It ends with
`ingest_status` set to `complete` or `failed`. Progress goes out on the SSE `curation`
frame with a `kind` of `ingest_started`, `ingest_progress` (with a `stage` of
`discovery` or `analyzing`), and then `ingest_complete` or
`ingest_failed`. Frames from a rescan also carry `phase: "rescan"`.

## Configuration fields

The experiment row keeps `name`, `description` and `path` (the beamtime folder,
set at create). `PATCH /api/experiments/:id` edits `name` and `description`; a
body with no patchable field, or with any other field (even alongside a name),
is a `400` that names the field.

Everything else is a column on the experiment's `setups` rows.

| Column | Set by | Editable later |
|--------|--------|----------------|
| `folder` | Create: the setup subfolder's name | No |
| `display_name` | `PATCH /api/setups/:id` (`name`) | Yes: shown instead of `folder`; an empty name or the folder name clears it |
| `included` | Create (the funnel's tick box), or `PATCH /api/setups/:id` | Yes: a hidden setup's samples leave the workspace and the scan skips it; showing a hidden setup starts a rescan, so frames that arrived while it was hidden come in. The last shown setup cannot be hidden |
| `mode`, `mode_source` | Create (the preview's detected mode, or the request's `mode`), or a scan that finds evidence | Not yet. `capillary`, `autosampler`, or `NULL` while undecided; a scan changes it only as described in [Setup modes](#setup-modes) |
| `root_dir`, `data_dir`, `analysis_dir` | Create | No |
| `image_pattern`, `metadata_pattern`, `integration_pattern` | Create | Not yet (see [Rescans](#rescans)) |
| `energy_kev`, `flight_path_m`, `beam_center_x`, `beam_center_y`, `pixel_size_um` | Create, or the scan's fill-only step | Not yet |
| `q_units` | Not yet | Reported as `A-1` with source `default` while `NULL` |

`GET /api/experiments` and `GET /api/experiments/:id` return the setups under
`setups`, ordered by `position`, with `name` set to the display name.

`PATCH /api/setups/:id` takes `name` (a string, or `null` to go back to the
folder name) and `included`. Any other field is a `400`. It logs a
`setup_updated` event and broadcasts it, and the palette's Settings page edits
both (undoable). In the workspace, samples of hidden setups are left out of the
list, Z, the progress count and every picker, and bridges or comparisons to them
stay stored but are not shown. A bridge is a merge of the same sample, so a
hidden mate still shares the call: a call edit on its visible mate changes it
too, and leaving that bridge needs the setup shown first. With more than one shown setup, the sample list
is grouped under a sticky header per setup.

### Migration from one experiment per setup

Before setups existed, each experiment row held one setup's configuration.
`open_db` migrates such a DB in place (`migrate_setups!`, `setups.jl`): each
experiment gets one setup copied from its old columns, a `path` equal to `data_dir` becomes the setup's
root folder, every load, sample and exposure gets that `setup_id`, the exposure
filename index becomes unique per setup, and the old configuration columns
(with `manifest_path`, `config`, `experiment_type`, `last_scan_tier` and
`consecutive_empty_ticks`) are dropped. Ids do not change. Merging several of
these experiments into one beamtime is a separate, one-time CLI step. The migration is one-way:
a build from before it cannot open the migrated DB, so copy `himalaya.db`
before the first start of an upgraded build.

## Setup modes

Each setup has a `mode`: `capillary` (a fixed-position capillary holder, the
default behaviour described in the rest of this document) or `autosampler`
(SSRL BL 4-2 solution SAXS with the autosampler). `NULL` means not decided yet.
`mode_source` records who chose it: `default` (no mode), `detected` (a scan or
create found it from the files), `user` (the request gave `mode`) or
`migration` (set by the upgrade). `GET /api/experiments[/:id]` returns both in
each setup.

**Detection.** `detect_setup_mode(data_dir)` (`setup_mode.jl`) lists `data_dir`
once. It answers `autosampler` when the folder holds at least one
`AT_SolSAXS_<label>_<scan>.txt`, a file of the autosampler whose trailing number
is a start scan, `capillary` when it holds at
least one `.tif` and no run file, and nothing when it holds neither. Image names
are not used: a capillary series named like `DOPE_pH7_37C_S0001_0_001.tif` can
look like an autosampler name and is detected as `capillary`.

**Where the mode comes from.** `GET /api/fs/preview` reports `mode` for each
setup (`null` for no evidence). `POST /api/experiments` and
`POST /api/experiments/:id/setups` take an optional `mode` (source `user`); an
unknown value is a `400`, and an absent one is detected from `data_dir` (source
`detected`). A detection with no evidence is stored as `NULL`. For an
autosampler setup, `analysis_dir` and any pattern the request leaves out are
filled from the autosampler layout. A setup found by a rescan
(`add_new_setups!`) takes its detected mode.

**Migration.** `migrate_setups!` adds the `mode` and `mode_source` columns if
the DB lacks them and, once (recorded as `setup_mode_v1`, after `setup_id` is
backfilled on exposures), sets `capillary` with source `migration` on every
setup that has exposures and no mode, in SQL, without reading the filesystem.
This also covers a DB from before the `setups` table existed. Setups without
exposures stay `NULL`. The columns come last in `setups`, in a fresh DB too, so
`SELECT *` returns the same fields in the same order everywhere.

**Layout.** `resolve_experiment_layout` returns the mode. For an autosampler
folder `analysis_dir` is always
`<root>/analysis/automatic_analysis/tot_files`, whether or not it exists yet (a
setup can be created after the run file but before the first integration), and
the patterns are `{name}<frame>.tif`, `{name}<frame>.prp` and `{name}_tot.dat`.
`<frame>` is the suffix of the lowest frame number among the folder's `.tif`
names ending `_<rep>_<frame>.tif`, for example `_0_001` or `_0_01` (beamtimes
differ in the frame number's width), and `_0_001` while there are no frames. The
beamline writes 10 to 20 one-second frames per acquisition, and `_tot.dat` is the
raw average of the acquisition's frames. Only the first frame is scanned, so an
exposure is one acquisition, its `filename` is the stem without the frame suffix
(`<label>_<well>_<S|B><scan>`), and `resolve_trace_path` finds
`<filename>_tot.dat` by the exact name. The scan ignores the other trace folders
(`sub_files`, `sastool`, `other_files`).

**Grouping.** An autosampler scan uses `group_autosampler`
(`autosampler_grouping.jl`) instead of `group_into_samples`. The capillary never
moves, so position would collapse every acquisition into a few samples.

- One sample per acquisition, holding that acquisition's first-frame exposure.
  The persist step never attaches an exposure to an existing sample in an
  autosampler setup (`reuse_samples` in `_scan_setup!`): an acquisition whose
  stem is not stored yet always gets a new sample. Scan numbers are not unique
  in every beamtime (the counter can restart), so two acquisitions can share one.
- Acquisitions with a well token are ordered by scan number. Samples are shot
  as a series: a run has one series label, and each well's plate coordinate
  (two digits and a letter, like a 96-well plate) is appended to it, so every
  well has its own name. A load is one run: a block of consecutive acquisitions
  with one series label. A new load starts when the series label changes, or
  when the acquisition's scan number is the start scan named by a
  `data/AT_SolSAXS_<label>_<scan>.txt` file. Without such files the label alone
  decides. One such file can be followed by several series before the next file
  appears (seen in a real folder), so the file does not mark one run.
- Stems without a well token (calibration, empty and test shots) all go in load
  `0`, never in a run's load.
- Indexes come from the files, not from counting, so a file that arrives late
  never shifts the others: a sample's `slot_index` is its acquisition's scan
  number (a stable hash of the stem when it has none) and a load's `load_index`
  is the scan number of its first acquisition. Loads are found or created by
  `load_index`.
- A sample is named `<label> <plate><well> (<S|B><scan>)`, for example
  `lipidA 01B (S0102)` or `water 01A (B0101)`. A stem without a well token keeps
  its stem as the name. This replaces the `(SNNPMM)` suffix of capillary samples.
- Rescans stay insert-only: new acquisitions add rows and never change existing
  ones. A run file that appears after its run's first acquisition was scanned
  splits the run from that scan number on: the earlier acquisitions stay in
  their load and the later ones get a new load. A metadata file that arrives
  late only delays its acquisition. If a run's first acquisition arrives after
  the second, the run is split the same way.
- Known limit: `load_index` is the first scan of a run, so if a folder's scan
  counter restarts, two runs can share a `load_index`. Ingest then reuses the
  first load row for both. Their samples stay separate, but that load's times
  and frame count are the first run's.

**No automatic analysis.** The scan stores no analysis for an autosampler setup:
no peaks, no `analysis_error` (so nothing in the problems list), no beam-off
deselection. `pending_analysis_ids` never reports its exposures, and
`ingest_has_work` still reports new images. The CLI `analyze` and the peak
model pass skip these setups, and `refresh_selection!` makes no average for
their samples. A sample's shown image is its one exposure, so the image strip
and the trace route work without analysis. The server refuses edits that would
analyze it: adding, excluding or removing a peak, adding or removing an
assignment member, a custom or speculative index, deleting an index,
`POST /api/exposures/:id/analyze` and `PUT /api/samples/:id/selection` all
answer 422 and write no event. Setting a call state is not refused: it cannot
start an analysis, and undoing a bridge between two autosampler samples
restores it. Labels, bridges within one mode, comparisons and drop or restore
still work.

**In the workspace.** The workspace picks its modules from the setup mode of
the open sample (`workspace/modules.ts`). A `capillary` or undecided setup has
all of them. An `autosampler` sample keeps the trace plot (zoom, plot view,
comparisons, subtraction), the image strip with Drop and Restore, labels and
bridges. It has no peak editing, no call keys or palette entries, no image
selection, and no review state. A module the mode lacks is not shown and the
rest move up: the tools column holds the Background section above the Image
section, open, in both the tall and the short layout. The column stays empty until the setups load, so the
first paint never shows modules the sample's mode then removes. The header's
progress block keeps its place, hidden, when no sample has a review state, so
it does not appear as data loads. A capillary sample and an autosampler sample cannot be bridged: the "Same sample as" page offers only samples of the open sample's mode (a setup with no mode counts as capillary), and `POST /api/samples/:id/bridge` answers 422 for a bridge, or a join of two groups, that would mix them. A mixed group stored before this rule is left as it is and stays listed. This arrangement is provisional: as modules are added to the autosampler mode, they go where they are sensible.

**Subtraction background (autosampler).** A sample's background defaults to the buffer at the start of its run: the buffer (stored exposure filename ending `_B<scan>`) with the lowest scan in the same load. A buffer, a run without a buffer and a capillary sample have no default. A buffer is a file whose stem has a well token and ends `_B<scan>`; files in load 0 (the catch-all for stems with no well) are never a buffer of a run and get no default. The default is computed on read (`default_backgrounds`, `subtraction_background.jl`) and never stored. A rescan is insert-only and a run's load is keyed on the scan of its first file, so a new buffer with a lower scan than the existing first buffer changes the default of the whole run only when it lands in the same load (its scan is after the run's first file); one scanned before the run's first file starts its own load and changes nothing. `sample_subtractions` holds only the user's choice: a row with a sample id, or a row with a null `ref_sample_id` meaning the user turned subtraction off. `GET /api/experiments/:id/compares` adds `backgrounds`, one entry per sample of an autosampler setup: `{sample_id, ref_sample_id, source, default_ref_sample_id, is_buffer, run_buffer_ids}`, where `ref_sample_id` is the effective background (null for none) and `source` is `default` or `user`. `subtractions` keeps listing only chosen samples. `PATCH /api/samples/:id/compares` takes `subtraction` as a sample id, `null` (none for an autosampler sample, cleared for any other) or `"default"` (delete the stored choice). For an autosampler sample the background is independent of its comparison chips: adding or removing a chip never clears it, and it uses no comparison slot. The workspace opens such a sample in Minus (the difference against the background) and shows the background in a Background module; the scale uses the same `compareScale` as a chip. A background in a hidden setup is treated as hidden and the plot stays Raw. `HimalayaDB.subtraction_background` returns only a chosen sample, and `nothing` for none or the default.

The sample list shows no unreviewed mark or call for it (an analysis-gap mark
such as "All dropped" still shows, without the word "Unreviewed"). It is not
counted in the header's "k of N reviewed" or on the home card
(`review_progress`), and `Z`, the landing sample and the prefetch of the `Z`
target skip it. The previous and next samples in the list are still prefetched,
whatever their mode. The home ribbon draws its span in a neutral shade and its
tooltip says nothing about review (the timeline's third value is `2`, against
`1` reviewed and `0` not). The selection controls follow the sample that owns
each image, so a bridged capillary mate keeps its selection controls beside an
autosampler sample, and the help lists the selection rows when any image in the
strip can be selected. When no
sample in the experiment has a review state, the header shows no progress and
the card shows no review count.

**Images the average kept.** An autosampler exposure is the beamline's average
of its acquisition's frames, minus the frames its software rejected, so it is not
a single image. The beamline's `sastool/<filename><frame>.log` (under
`<root>/analysis/automatic_analysis`, `<frame>` from the setup's image pattern)
holds the counts. `exposures.frames_total` and `frames_rejected` (the last two
columns, NULL when unknown) take them from the `Sample frames:N Rejected:M`
line. The log of a run's first buffer has no summary, so its `frame #` lines are
counted instead, once the log has been unchanged for a minute; that count is not
revised afterwards. Only a log with a `Sample frames:` line, or with no summary
lines at all, gives a count. A log with no readable counts, or no log
(calibration shots), leaves both NULL. Capillary
setups are never read, and a stored value is never overwritten. `backfill_frame_counts!`
fills unset rows at the end of every `scan_and_group!`, on every
`POST /api/experiments/:id/scan` and once at boot, so rows ingested earlier
are filled without new files and a log that appears later is picked up on the
next scan. The Image panel then reads `9 of 10 images, beamline average` (or
`10 images, beamline average` when none were rejected) and an exposure of `9 × 1.0 s`; with
NULL counts it reads `Single` and `1.0 s` as before.

**Settling a setup.** At the start of every `_scan_setup!`:

- A setup with mode `NULL` or mode `capillary` with source `detected` is
  promoted to `autosampler` (source `detected`) when the folder now has a run
  file, no stored exposure filename carries a well token, and no stored exposure
  is reviewed (`reviewed_at` set). Promotion rewrites `analysis_dir` and the
  three patterns to the autosampler layout and renames the stored exposure whose filename ends with the layout's first-frame
  suffix (a calibration shot taken before the first run file) to its acquisition
  stem, so the layout does not ingest it a second time. Its other frames keep
  their names.
- Otherwise a setup with mode `NULL` and a detected `capillary` stores
  `capillary` with source `detected`. With no evidence it stays `NULL`, scans as
  capillary (there is nothing to scan) and is detected again on the next scan.
- A setup with source `user` or `migration`, or a `capillary` setup with
  exposures that carry well tokens or are reviewed, is never changed by a scan.
- An `autosampler` setup that has no exposures yet gets its layout derived
  again on every scan, which covers a run file that exists before any frame does
  and a frame width that is only known once frames arrive.

A hidden setup (`included = 0`) is not scanned, so it settles when it is shown.

## File patterns

Each pattern contains `{name}` exactly once; it stands for an exposure's stem.
A pattern without it would point every frame at the same file, so it is
refused everywhere one is accepted: `400` from `POST /api/experiments`,
`POST /api/experiments/:id/setups` and `GET /api/fs/manifest`, an
`ArgumentError` from `create_setup!`, and an inline error in the funnel's
Details that blocks adding a ticked setup (an unticked one is sent without
patterns).

| Pattern | Default when the column is `NULL` | Matched in |
|---------|-----------------------------------|------------|
| `image_pattern` | `{name}.tif` | `data_dir` |
| `metadata_pattern` | `{name}.prp` | `data_dir` |
| `integration_pattern` | `{name}.dat` | `analysis_dir` |

The scan enumerates the images, takes each stem, and builds the metadata and
integration paths by putting the stem into those patterns. An exposure exists
for every image, whether or not its metadata or integration file is present.

For a capillary folder the resolver recognises two layouts by probing one image stem (an autosampler folder uses the layout in [Setup modes](#setup-modes) instead):

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
No view draws them yet.

## Scan and grouping

`scan_and_group!(db, experiment_id)` is the only code that inserts exposures. It
runs once per included setup (`_scan_setup!`), reading `data_dir`,
`analysis_dir` and the patterns from the setup row, then:

1. **Discovers** exposures with `scan_directory` (`grouping.jl`): one per image
   in `data_dir`, with its metadata file parsed (`parse_prp`, `prp.jl`) and its
   integration file located. If no file matches the image pattern, it returns
   `status = :empty` and writes nothing. When the folder has any metadata files,
   an image whose metadata file has not appeared yet is left for a later pass
   instead of being grouped without a position or time. A folder with no
   metadata files at all is still ingested.
2. **Fills geometry** that is still unset (see above).
3. **Groups** exposures with `group_into_samples` (capillary setups; see [Setup modes](#setup-modes) for autosampler setups): it splits them into loads at
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
   `I_2`, the beamstop photodiode, `incident_monitor` from `I_0`, `counting_time`).
   `I_signal` is not read: it copies whichever counter the `.prp` names, usually
   the APD, which reads close to the incident beam. Each setup's trace
   normalisation comes from the nearest `integ.mpp` at or above its analysis folder,
   searching no higher than the experiment's own folder
   (`setups.trace_normalisation`: `i2`, `i1`, `none` or `unknown`). It is re-read with the
   monitors on every scan and after new setups are added, and a changed `integ.mpp`
   updates the stored value. An unreachable analysis folder keeps the stored value.
   SSRL writes it at `<setup>/analysis/integ.mpp`, two levels above `tot_files`, and
   lines starting with `#` are skipped. Only autosampler setups use it for the buffer
   scale (scale 1 on `i2`, the `I_2` flux ratio on `none`); capillary setups match at high q.
   `backfill_monitors!` reads them from the stored `prp_path` of every exposure
   still marked `monitors_read = 0`, on every scan request and at the start of
   every pass, so a finished beamtime ingested before the monitors
   were parsed picks them up with no new files. Each file is read once, and a
   stored value is never overwritten, with one exception: the
   `exposures_transmitted_from_i2_v1` migration clears every stored
   `transmitted_counts` once, because it held `I_signal`, so the backfill reads
   `I_2` in its place.
5. **Analyzes** (capillary setups only), outside the write transaction, every exposure that is new,
   never analyzed, whose integration file has appeared since it was last
   tried, or whose stored `analysis_version` differs from `ANALYSIS_VERSION`
   (`pipeline.jl`: the core Himalaya version plus a local counter, bumped when
   the analysis changes). An analysis from before versions were recorded has
   no `analysis_version` and is re-analyzed too: the migration that adds the
   column leaves it empty rather than claiming the current version. A reviewed
   exposure (one a user has curated, assigned, called or accepted) is never
   re-analyzed for a version change, only when its trace file changes. At SSRL BL 4-2 the image and metadata file are written
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
(for example "2 images waiting for a trace · 1 analysis failed") and its
Details button opens the list; picking a row opens that sample. A frame
waiting for its trace says "Waiting for this image's trace." in the plot, not
"No integrated trace", and the plot refetches the trace when the next pass
completes.

It also returns `setups`: each folder under the beamtime that would be a new
setup but whose file layout `resolve_experiment_layout` does not recognise
(`unrecognized_setup_roots`, `{root, folder}`). A rescan skips these, so the
notice counts them ("1 setup folder not recognised") and the list names them,
with an "Add setup folders" link that opens New experiment on the beamtime
folder. There the funnel shows the beamtime as already added, keeps its
setups locked, and lets the user tick the new folder, type its patterns in
Details, and add it with `POST /api/experiments/:id/setups`.

`POST /api/experiments/:id/setups` takes `{setups: [...]}` in the same shape and
with the same checks as the beamtime create (each root inside the experiment's
`path`, each data folder unused, at least one ticked). It adds them as the next
positions under the write lock, starts a rescan and answers
`202 {status: "scanning", experiment_id, setup_ids}`. An experiment whose `path`
is its own single setup folder cannot take more setups (`400`).

## Rescans

Every pass goes through `run_ingest!` (`server.jl`), which holds a per-experiment
lock, so one experiment never has two passes at once; a second caller waits and
then runs its own, usually empty, pass. It sets `ingest_status` to `scanning`,
back-fills beam monitors, adds new setup folders (`add_new_setups!`,
`beamtime.jl`), runs `scan_and_group!`, suggests role tags, and ends
`complete` (stamping `last_scanned_at`) or `failed`. Every pass after the first
sends `phase: "rescan"`.

Role tags come from `suggest_labels!` (`label_suggestions.jl`, word lists in
`ROLE_WORDS` and `SOLVENT_WORDS`). The name is split into words, each also tried
without trailing digits:

- `agbe`, `agbeh`, `behenate`, `calibrant` or `lab6` suggests Calibrant;
- `empty`, `emptycap` or `emptycapillary` suggests Empty capillary;
- `buffer`, `buf` or `blank` suggests Buffer;
- `hepes`, `pbs` or `tris` suggests Buffer only when every other word is a
  code (`S01`, `2`, or the `S01P02` slot code ingest adds) or another buffer word, so `HEPES_S01` is a buffer and
  `MO_HEPES_S01` or `DOPC_in_PBS` is not.

The tag is *suggested*, with source `ingest`. A sample that already has any role
tag, suggested, confirmed or removed, is left alone, so a rescan never brings
back a suggestion the user removed or changed. The step runs at the end of every
pass, after `scan_and_group!`. If it fails, the pass still ends `complete` and
the failure is only logged. Experiments ingested before suggestions existed get
them once, when the database is opened (`label_suggestions_v1`). A suggestion
never picks a subtraction buffer by itself. It only puts tagged buffers first in
the compare picker, where confirmed buffers come before suggested ones. Ranking
buffers by their traces (no peaks, a flat trace, a low floor) is not built yet
(#368).

The same step suggests a flight-path label (`suggest_flight_path!`). Each
exposure stores the `.prp` pipe length, which is the flight path, as
`exposures.pipe_length_m` (read at scan time, or by the beam-monitor back-fill).
When every exposure of a sample that has a pipe length rounds to the same
tenth of a metre, the sample gets the suggested label `<n.n> m`, for example
`1.7 m`, with source `ingest`. Frames that disagree, or no pipe length at all,
suggest nothing. A sample that already has any `<number> m` label (lowercase
m, matched case-sensitively, so a molar `1 M` doesn't count), in any state, is
left alone, so removing or renaming it sticks. Databases from before
this column re-read their `.prp` files once (`exposures_pipe_length_v1` resets
`monitors_read`), and the start-up back-fill then suggests the labels.

A new setup folder is any folder `setup_roots` finds under the experiment's
`path` that is not one of its setups yet and whose data folder no other
experiment uses. It waits until `resolve_experiment_layout` detects its layout
(an integrated trace, or the SSRL `analysis/automatic_analysis/tot_files`
folder), so a folder whose images arrive before anything else is not stored
with guessed patterns (until then it is listed under [Problems](#problems)).
It is then added as the next `position`, ticked, with
that layout; the scan then fills its geometry. Setups
the user unticked stay setups, so they are never re-added. An experiment whose
`path` is one of its own setup folders (a single setup was picked, or a DB
migrated from one experiment per setup) never adds folders. An unreachable
`path` adds nothing.

`ingest_has_work(db, id)` (`ingest.jl`) says whether a pass would change
anything: a new setup folder, an image (with its metadata file, when the folder has any) that has no
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

## CLI reference

The CLI does not create experiments. It works on experiments already in the
central DB (`HIMALAYA_DB_PATH`, else `~/.himalaya/himalaya.db`).

```
himalaya analyze [-e <experiment> | --all] [-s <sample>]
    Run peak-finding and indexing for every frame that is not dropped.

himalaya show [-e <experiment>] -s <sample>
    Print the stored peaks, indices and active group for one sample.

himalaya serve [--port N] [--host H]
    Start the web server. Creates an empty DB if none exists, so a new user
    can create the first experiment through the funnel.

himalaya merge <experiment> <experiment>... [--name N] [--path P] [--apply]
    Merge experiments made one per setup into one beamtime. Without --apply
    it prints the plan and writes nothing.
```

`merge` keeps the first experiment named and moves the others' setups, loads,
samples and frames into it, then deletes the emptied experiments. No other id
changes, so phase calls, labels, bridges and the event log carry over. The
result's `path` is the folder the setups share (or `--path`) and its name
replaces the kept experiment's name: that folder's name, or `--name`. Setups
are ordered by folder name. It refuses fewer than two experiments, one listed
twice, one that is scanning or analyzing, two setups with the same folder name,
setups that are not all in one folder (unless `--path` names it), and a `path`
another experiment already has. Stop the server before applying it, so no
open page holds stale experiment ids.

`-e` takes an id, a name or a path. A path matches the experiment's folder whether
it is typed through a symlink or with a trailing slash. When it is left out, `analyze` and `show`
use the only registered experiment and fail with a list if there are several.

## Further reading

- [packages/HimalayaUI/src/routes_fs.jl](../packages/HimalayaUI/src/routes_fs.jl): the funnel probes
- [packages/HimalayaUI/src/routes_experiments.jl](../packages/HimalayaUI/src/routes_experiments.jl): create, `PATCH`, scan
- [packages/HimalayaUI/src/ingest.jl](../packages/HimalayaUI/src/ingest.jl): `scan_and_group!`, `ingest_has_work`, `pending_analysis_ids`
- [packages/HimalayaUI/src/grouping.jl](../packages/HimalayaUI/src/grouping.jl): discovery, loads, slots, sample names
- [packages/HimalayaUI/src/geometry.jl](../packages/HimalayaUI/src/geometry.jl): geometry derivation
- [packages/HimalayaUI/src/config.jl](../packages/HimalayaUI/src/config.jl): `integration_pattern_for`
