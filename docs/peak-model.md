# Per-database peak model

Status: implemented (all four delivery steps). Packaging, provenance and pinning, the training service and model-written auto peaks are built: the features, labels,
model and cross-validation live in `packages/HimalayaPeakModel`, the database
loader is `HimalayaUI.labeled_exposures` (reading review snapshots), `tools/peak-classifier/` (#486)
is the command-line tool on top of both, and the server trains and records runs in the background
(see "Training service (step 3)"), promotes a passing run and writes its peaks (see "Model-written peaks (step 4)"). This document describes how HimalayaUI trains such a model for its own
database, keeps it up to date as users curate, and uses it.

## Goal

Each HimalayaUI database learns from its own curations which candidate peaks
are real, without data or models leaving the machine. Training, evaluation and
updates run in the background with no user action, and the user experience
does not change: users see better auto peaks, not a model. A database with too few
curations, or a model that does not beat the hybrid detector on that
database's own held-out data, falls back to the hybrid.

## Constraints from the current backend

- `serve` serves one database, so "per database" means per server process.
- Background work: `start_background!` spawns one task at boot
  (`backfill_all_monitors!`, then `resume_ingests!`). The only periodic work is
  the idempotency GC timer. There is no job queue. Writes serialise on
  `_DB_WRITE_LOCK`; GET requests read from pooled read-only connections.
- Analysis memoisation: `ANALYSIS_VERSION` is a constant. An exposure whose
  stored `analysis_version` differs is re-analysed by the next ingest pass,
  and a curation route on such an exposure runs a full `findpeaks` inside the
  request until that happens.
- Curations are stored by q and re-attached to auto peaks by
  `peak_matching.jl`, so they survive a change of peak set. Auto peak ids are
  kept across re-detection when q matches (`diff_update_auto_peaks!`).
- There is no settings or key-value table, and the schema avoids BLOBs. Files
  next to the database have a precedent (`thumb_cache_dir`), and `open_db`
  makes the database group-writable for shared deployments.
- The tool depends on HimalayaUI by path, so HimalayaUI cannot depend on the
  tool as it stands.

## The model's role

The active model writes `auto_peaks` for exposures nobody has reviewed yet.
`analyze_exposure!` runs the hybrid detector as before, then the model over the
candidate pool, and stores the model's peaks. Users see the same workspace and
the same peak marks, and curate exactly as today: add peaks, remove added
peaks, exclude and un-exclude auto peaks. Curations are stored and matched the
same way whichever detector produced the auto peaks. Only which auto peaks
appear changes, and those same curations are the labels the next training run
learns from.

Measured basis: see "Evidence". The model ties the hybrid on phase calls and
finds markedly more of the peaks users otherwise add by hand, so the gain is
less manual peak work.

## Reviewed exposures are pinned

An exposure is **pinned** (also called reviewed) once one of these events
touches it: a peak add, remove, exclude or un-exclude; an assignment add,
remove or state change; a speculative index creation; a sample call (which pins every exposure of that
sample, including calls made by a bridge); or an exposure accept or reject. The
retired index confirm and un-confirm kinds are never emitted and are not in the
list. The event sets
`exposures.reviewed_at` in the same transaction, so the flag cannot lag the
event. Pinning is broad on purpose: it only stops re-detection. An exposure is
**labelled** only when the user acted on its own peaks or indexing, which is
narrower (below). The training trigger's label kinds are the labelling events,
not every pinning event.

A model promotion never re-detects a reviewed exposure. `persist_analysis!`
re-attaches assignment members by phase and basis, so a changed peak set can
drop a user's assignment, and open clients would get the assignment-dropped
toast. Pinning prevents both. Curation routes never re-detect because of the
model either: the curation itself pins the exposure with the peaks the user
saw.

An `ANALYSIS_VERSION` bump does not re-detect reviewed exposures either: the
ingest pass skips a reviewed exposure whose version is stale, and
`analyze_exposure!` honours a reviewed exposure's stored hashes whatever its
version, so a curation never triggers re-detection. Only a changed trace file
re-detects a reviewed exposure. Re-detecting would change the peaks behind the
user's curations and could drop assignment members (`persist_analysis!`
re-attaches them by phase and basis). Unreviewed exposures are re-analysed on a
bump as before. Labels are unaffected either way, because training reads
review snapshots, not live auto peaks (below).

## Provenance bookkeeping

Each auto peak records where it came from, set by `diff_update_auto_peaks!`:

- `auto_peaks.found_by`: `hybrid`, `model` or `both`;
- `auto_peaks.model_id` (NULL when no model was active).

Each snapshot also records whether it is **direct**: 1 when the event targets
the exposure itself and is a peak add, remove, exclude or un-exclude or an
assignment add, remove or state change; 0 when it was written because of a
sample call or bridge on other frames, or an accept or reject. Training labels
only exposures with at least one direct snapshot (plus form-factor ones). A
sample call marks many frames the user never looked at, and treating all their
peaks as accepted gave far more positives and weak negatives than the user
curated. A bridge-induced call on another sample pins and snapshots those
frames but never makes them labelled.

Labels do not come from live `auto_peaks`, which a release can re-detect. Each
reviewing event instead writes a **review snapshot** in the same transaction:
the exposure's auto peaks at that moment (q, `found_by`, `model_id`), the
detector on screen (`ANALYSIS_VERSION` plus model id), and the curation.
Training reads only snapshots, using each exposure's latest **direct**
snapshot for the peaks the user saw (a later indirect snapshot, such as a
sample call after a re-detection, is ignored) and all its curations as they
stand now.

Labels for training (all by q). They use the app's tolerance function,
`peak_match_tol`, with the grid step taken from the trace and an
any-candidate-within-tolerance test, whereas the app attaches each curation to
its single nearest auto peak using `exposures.q_step`. Rules are applied in
order, first match wins:

| Evidence | Label |
|---|---|
| nearest candidate to an `add` curation | positive |
| nearest candidate to a snapshot peak left alone, `found_by` = `hybrid` or `both` | positive |
| other candidates within tolerance of a positive | dropped |
| candidates within tolerance of a snapshot peak left alone, `found_by` = `model` | dropped |
| candidates within tolerance of an `exclude` curation, whether or not an auto peak still matches it | strong negative |
| remaining candidates on a form-factor exposure | strong negative |
| remaining candidates on a labelled exposure | weak negative |

Kept peaks count as positives on every labelled exposure, including ones with
no add or exclude (for example only an assignment edit). Today the tool counts
them only on exposures with at least one add or exclude, so this widens the
positive set.

Excludes are labelled at their own q (every candidate within the match
tolerance of an exclude is a strong negative): after a re-analysis many excludes match
no auto peak (the current detector does not find those peaks), and dropping
them would lose real negatives. A model-only peak the user left alone is
dropped rather than left to the weak-negative rule, which would otherwise
make it, and its neighbours, a negative. Exposures reviewed before snapshots
existed get one snapshot from their current auto peaks when the migration
runs, all `found_by = hybrid`, with `direct = 1` when the exposure has a
peak curation or a logged event of a direct kind on it (including the legacy
`index_confirmed` and `index_unconfirmed`, which mark pre-Plan-A human
confirmations) and 0 when it was pinned only by a sample call, a status, or an
assignment member with no such event. Members alone can be the auto guess that
`migrate_assignments!` seeded, so they do not count.

Implementation (step 2). The pin and snapshot are written by
`record_event_reviews!`, called from `apply_event!` after the dispatcher (not
from `update_view_for_event!`). `exposures.reviewed_at`, `peak_reviews` and
`auto_peaks.found_by`/`model_id` are bookkeeping, not views: the log cannot
reproduce them, so `rebuild_views_from_log!` leaves them alone (see
`docs/event-log.md`). `peak_reviews` holds `exposure_id`, `event_id` (NULL for
backfilled rows), `created_at`, `shown_detector` (`ANALYSIS_VERSION`),
`shown_model_id` (NULL until models exist), `backfilled`, `direct` and
`peaks_json`.
The backfill pins exposures that have a curation, an assignment member, an
accept or reject status, or a sample call (all frames of the sample), with
the first matching event time or the migration time. `labeled_exposures`
includes exposures with a direct snapshot and reads each one's latest direct
snapshot. Kept peaks with `found_by = model` go in `LabeledExposure.kept_model`,
and `label_candidates` drops them and their neighbours. `LabeledExposure.labelled`
is true for exposures with a direct snapshot.

Bridge call changes emit `sample_call_set_state`, and a `keep_call` bridge
also emits `assignment_set_state` on the losing side's frames. Both pin and
snapshot those frames; the bridge passes `review_direct = false`, so none of
them becomes labelled. The share-bridged-labels and share-bridged-calls
migrations run with pinning off: they are system changes, and the backfill that
runs after them pins the resulting calls as indirect. The provenance columns
are added before them (`migrate_peak_provenance_schema!`).

Growth: every reviewing event appends a full `peaks_json`, so a sample call
adds one row per frame of the sample, and bridged mates add more. Nothing
prunes `peak_reviews`. If that matters, store a new snapshot only when its peak
set differs from the latest one for the exposure.

`diff_update_auto_peaks!` keeps a re-matched row's `found_by` and `model_id`
unless the caller passes `found_by` (then both are written, `model_id` NULL
when omitted); inserts default to `hybrid` and NULL.

## Training

- **In a child process.** DecisionTree's `build_forest` and the feature
  extraction use `Threads.@threads` over every default thread, and their
  allocations trigger stop-the-world collections that would stall requests.
  Training therefore runs in a separate `julia` process, started with the
  server's sysimage so it does not spend minutes compiling, that opens the
  database read-only (`mode=ro`), trains and cross-validates, and prints the
  fitted model, recipe and metrics as JSON. The server writes the result under
  `_DB_WRITE_LOCK`.
- **Trigger.** A watermark: the largest `user_actions.id` among the label
  kinds (the direct labelling events above). Training is due when at
  least N new label events have arrived since the last run and no label event
  has arrived for T minutes. Proposed N = 20, T = 10.
- **Single flight.** One training process per server. A trigger during a run
  marks the next run due instead of starting another.
- **Feature cache.** Features depend only on the trace, the feature version
  and the core version, so they are cached per trace hash as raw Float32
  arrays in files next to the database (derived data, safe to delete; a raw
  array runs no code when read). Whoever can write that directory can alter
  training inputs, which is the same trust as write access to the database. For
  a hypothetical database of 10,000 exposures with 30 candidates of 28
  features each, that is about 34 MB as Float32, against well over 100 MB as
  JSON inside the database.
- **Cold start.** No model is trained until the database has at least 100
  add and exclude curations across at least 30 curated exposures, in at least
  3 experiments (or 5 samples, grouping by sample, when it has fewer
  experiments). Kept hybrid peaks are plentiful and say little, so they do not
  count towards this. Until then the hybrid alone is used, and the run's row
  records why.

## Training service (step 3)

What is built, and where it differs from the sections above.

**Child process.** `himalaya train-peak-model [--db PATH] [--trees N] [--no-cache]` opens the
database read-only (`HimalayaUI.open_readonly_db`, also used by the tool), runs
`run_peak_training` and prints one JSON object on stdout (status, metrics, recipe, fitted, watermark,
counts); progress goes to stderr. The server starts it with `Base.julia_cmd()`, which carries the
sysimage the server runs with, and the active project. The child gets at most
`min(server default threads, CPUs / 2)` threads and, on Linux when `nice` exists, runs under `nice -n 10`, so a run does not take
every core from requests and ingest. Its stdout and stderr go to
temporary files, so a crash leaves the end of stderr (about 1700 characters) in a `failed` row.
`PEAK_TRAINING_COMMAND` replaces the command in tests. An `atexit` hook and `stop_peak_training!` kill a running child. `atexit` does not run when
systemd stops the server with SIGTERM, so the unit must keep `KillMode=control-group` (the systemd default,
set explicitly in `deploy/himalaya.service`), which also stops the child.

**Trigger.** The watermark is the largest `user_actions.id` among `DIRECT_REVIEW_KINDS` on exposures (the
labelling events; `set_exposure_status` and sample calls do not count, unlike the review-snapshot
kinds). A run is due when at least 20 such events are newer than the largest watermark in
`peak_models` (every status counts, so a crashing run waits for new labels) and the newest has been
quiet for 10 minutes (`PEAK_TRAINING_MIN_EVENTS`, `PEAK_TRAINING_QUIET_SECONDS`). A timer checks
every 5 minutes, first after 60 seconds, and runs only when the server resumes work
(`start_background!(...; train = resume)`), so tests that boot with `resume = false` start nothing.
The row's watermark is the child's, read at its start, so labels that arrive during a run stay due.
One child at a time: a due check during a run returns `:busy`. `stop_peak_training!` sets a `stopped` flag
(cleared only by `start_peak_training_timer!`), under which a check that was already in flight returns `:stopped` and spawns
nothing; the child is registered under the same lock that stop takes, so stop kills it instead of waiting for it. The watermark
a run started from is also kept in memory, so if writing its row fails the next check does not rerun it until new label events arrive. There is no extra flag, because dueness is read from
the table, and the finished run's watermark leaves the newer events due. The server writes the row under
`_DB_WRITE_LOCK` on its own connection. No SSE frame, no `user_actions` row, no UI.
An in-memory database never trains (the child needs a file).

**Cold start.** `cold_start_check` counts add and exclude rows of `peak_curations` on non-rejected file
exposures. It needs 100 curations on 30 exposures, and 3 experiments, or else 5 samples (then the
folds group by sample). Otherwise the child trains nothing and returns `rejected` with one reason per shortfall.

**Dataset.** The run loads every non-rejected exposure with a readable trace (`include_unlabelled`), not only
the labelled ones, so the peak caps below are measured over "all exposures with a trace". Unlabelled
exposures carry no labels; they are never trained on. `oof_proba` predicts a group that has
no labelled candidates (an experiment or sample nobody curated) with one model trained on everything, which
is what leaving it out would give, instead of fitting one forest per group.

**Feature cache.** `HimalayaPeakModel.DiskFeatureCache`, in `<database dir>/cache/peak-features/<database file stem>/` (so databases in one directory do not prune each other's files),
one file per exposure named `<sha256 of the q, I, σ arrays>-<first 12 hex of sha1(feature version)>.f32`
(a hash of the loaded arrays rather than `exposures.trace_hash`, which can lag the file). Layout:
8-byte magic, Int32 `n` candidates, `p` features and `h` hybrid peaks, then `n·p` Float32 (column-major), `n`
Int32 candidate indices and `h` Int32 hybrid indices. The hybrid indices are stored because the gate needs them.
A short, malformed or wrong-width file is a miss and is rewritten. A changed trace or feature version is a new
file, and each run deletes every `.f32` file it did not use and `.tmp` files over an hour old (left by an interrupted write). An unwritable or racing
directory is ignored per file (writes and deletes). The training run always rounds features to Float32, with the cache or with `--no-cache`, so every run builds the same dataset.

**Gate.** All figures are out of fold. The threshold is chosen on the same out-of-fold predictions and label set the checks then score, so
add recall is the best of the 13 thresholds (optimistic), and the two peak-cap checks pass by construction whenever some threshold fits the caps; the excludes and
phase checks are not affected by the choice. `peak_summary` scores a detector's peaks, and `evaluate_gate` applies the checks to
the model's and the hybrid's summaries (a check with no data is skipped and passes).

- Peaks per exposure ≤ hybrid + 0.5 over all loaded exposures, and ≤ hybrid + 0.1 over form-factor ones (constants `GATE_PEAKS_SLACK`, `GATE_FF_SLACK`).
- Phase calls. The exposures scored are labelled ones (a direct snapshot with no model shown) whose assignment is `indexed`
  (and sample call `indexed`, the default), with the call being the `(phase, basis)` of the members. Members on `speculative`
  indices (the custom indices users build) always count, being user-made. Members on `auto` indices (the auto-group picks, which are fitted from a
  detector's peaks) count only when the exposure has a user `assignment_add`, `assignment_remove` or `assignment_set_state` event
  (`phase_auto = :user_events`, the default; `:all` keeps every auto one). The stored calls do not record confirmation beyond that, so an
  `auto` call with such an event is the closest signal to a confirmed one. On a local real database the event-limited set kept most of the
  exposures and did not change the gate outcome, so the limit is used. A pick is correct
  when the first `auto_group` pick over the detector's peaks, with excludes removed, adds not applied and
  sharpness from `sharpness` with the pipeline's fallback (median positive), has the phase and a basis within 2% of one of the calls.
  Model correct count must be at least the hybrid's.
  Only reachable calls are scored (`phase_targets = :reachable`, the default; `:all` scores every call). A call is reachable when `indexpeaks`
  yields a candidate of the called phase with a basis within 2% of the call, over the user's own final peaks (kept snapshot peaks and
  adds snapped to the grid, excludes removed, the same sharpness rule) or over the peaks of either scored detector (the model at the chosen
  threshold and the hybrid, excludes removed). Calls no scored peak set can produce, for example a custom index anchored on fewer peaks than
  indexing needs, only dilute the metric, so an exposure with no reachable call leaves the targets. An exposure with at least one keeps its full
  call list, so its phase key is the same in both modes and a hit on any call counts. A call whose phase exists but is out-ranked stays in.
  The model and the hybrid are scored on identical targets. The metrics JSON `comparison` holds `phase_targets` (the scored count; before
  this change it counted every confirmed call, so recorded runs before and after are not comparable), `phase_targets_mode` and
  `phase_target_counts` (`all`, `reachable` and `per_phase` counts).
- Add recall (hits over adds, adds in the top 7% of q excluded and counted as `adds_in_band`) must be at least the hybrid's, and the share
  of excludes still detected at most the hybrid's, on exposures whose latest direct snapshot had `shown_model_id` NULL. Today that is every labelled exposure.
- The threshold is chosen from 0.30 to 0.90 in steps of 0.05 as the one with the highest add recall within the peaks caps (higher on ties).
  When none is within the caps the run is `rejected` with that reason.

A passing run trains the final model on all rows. A rejected run trains nothing.

**Metrics JSON.** `cv` (the `cross_validate` tables per threshold, hybrid and candidate-pool rows, AUC), `comparison` (chosen threshold, model and
hybrid summaries, phase target count), `sweep` (per threshold: peaks, form-factor peaks, adds, within caps), `gate` (passed, checks, reasons, caps) and `seconds`.

Promotion, the 25% label-growth rule, `shown_model_id` and prediction in `analyze_exposure!` are step 4.

## Swapping the model

The model is a recipe (`HimalayaPeakModel.Recipe`) made of three parts that
can each change without touching the rest:

- **Learner** (`Learner`/`Fitted`: `fit`, `predict_proba`, `to_dict`,
  `from_dict`). The package has `RandomForestLearner` and `LogisticLearner`.
  A new learner is a new type plus one `register_learner!` entry (from the
  defining package's `__init__`; a duplicate kind is an error). Nothing else
  changes: training, cross-validation, the gate and storage take any learner,
  and old rows keep loading by their `kind`.
- **Feature set** (`FeatureSet`: `feature_names`, `feature_version`,
  `candidates_and_features`). `DefaultFeatures` is today's. A new feature set
  has a new feature version, so cached features are rebuilt and models trained
  on the old one are retrained from the stored labels.
- **Recipe**: feature set, learner, label options (`weak_keep`), threshold and
  NMS window, with a stable `recipe_id` (hash of its dict). Changing any
  setting is a new recipe and a new id.

Scope: a recipe describes a per-candidate classifier (candidate rows, a
probability threshold, NMS). A model of another shape, such as a whole-trace
model, would need `Recipe` and `TrainedModel` to become abstract with their own
`train` and `predict_peaks`.

A swap does not decide anything by itself. A model trained from a new recipe
is a pending candidate and the promotion gate below decides whether it replaces
the active one, using the same hybrid comparison for every recipe.

## Evaluation and promotion

A label is biased against the detector that was on screen when it was made:
an add exists because that detector missed the peak, and an exclude because it
showed one. While only the hybrid has been shown, every label is biased
against the hybrid, so add recall and excludes-detected favour any candidate
by construction. The binding checks are therefore the ones that do not depend
on who was shown:

- **Peaks per exposure**, over all exposures with a trace, must stay within
  the hybrid's plus 0.5, and peaks per form-factor exposure within the
  hybrid's plus 0.1. These cap how many extra peaks a model may add for its
  higher recall. The values are starting points to set from the first local
  runs.
- **Phase calls**: top-1 calls (first pick of `auto_group`) on exposures with
  a confirmed call, with excludes applied and adds not applied, must not fall
  below the hybrid's. Applying adds would hand both detectors the user's
  final peak set and hide any difference.
- **Add recall and excludes-detected** are reported, with each label's shown
  detector, so labels made under the current hybrid, older detectors and
  models can be told apart. They must not be worse than the hybrid's on labels
  made while the hybrid or an older detector was shown.

Every comparison is against the hybrid, scored out of fold (leave one
experiment out). There is no model-against-model test: two runs of the same
recipe differ only in data, threshold and seed, which the hybrid comparison
already covers. A new model replaces an active one only when labels have grown
by at least 25% in adds and excludes since the active model's watermark and
the new one passes the hybrid comparison, so models do not churn.

The threshold is the one that maximises add recall within the peaks caps. Adds
in the high-q band the detector trims (the top 7% of the q range for the
hybrid, 5% for the persistence default) are reported
separately and left out of add recall. A run that fails is recorded and
production stays.

## Storage

- Table `peak_models` (migration `migrate_peak_models!`, idempotent): id, created_at, status (`passed`,
  `rejected`, `active`, `retired`, `failed`), `label_watermark`, `counts_json`, `threshold` and
  non-maximum-suppression window `nms`, `metrics_json` (the cross-validated tables, the
  hybrid comparison and the gate outcome with its reasons),
  `feature_version`, `core_version`, `recipe_json` (the recipe's `to_dict`), `recipe_id`,
  and `fitted_json`. A run that passes the gate is `passed` (the design's earlier `pending`);
  `failed` is a crashed or unreadable run, whose `metrics_json` holds `{"error": <stderr tail>}`.
  Step 3 writes only `passed`, `rejected` and `failed`; promotion (step 4) moves rows to `active` and `retired`.
- The fitted model is `to_dict` of the learner's `Fitted` object, stored as JSON
  in a TEXT column, only for the newest `passed` row and any `active` row; it is NULL
  on `rejected`, `failed` and `retired` rows, and inserting a `passed` row clears it from older `passed` ones. For the random forest each tree is nested split
  nodes (exact Float64 thresholds) and leaves holding the majority class and the
  training count per class. JSON loads without running code (a `Serialization`
  file runs code when loaded, and the directory next to the database is
  group-writable) and does not depend on the Julia version. Loading dispatches
  on the dictionary's `kind` through a registry and rebuilds DecisionTree 0.12's
  own tree structs, so prediction is DecisionTree's code. That ties the stored
  forest format to DecisionTree 0.12's struct layout; `[compat]` holds the
  bound, and the round-trip tests (identical forest and per-tree
  probabilities) must pass on any upgrade. The recipe row is enough to retrain
  the same model from the stored labels.
- The feature version includes `pkgversion(Himalaya)`, because candidates and
  the detector-agreement features come from `findpeaks`. After a core upgrade
  the server retrains from the stored labels at boot and runs the gate, so a
  release re-analyses unpinned exposures once, not model → hybrid → model.
- Rebuild migrations that list `auto_peaks` columns explicitly must carry
  `found_by` and `model_id`.

## Analysis with a model

- `exposures.peak_model_version` records which model produced an exposure's
  peaks. `analysis_version` stays equal to `ANALYSIS_VERSION`, so the ingest
  queries, `ingest_has_work` and the experiment-open scan never see a
  promotion: no "Analyzing" badge and no ingest notice.
- One function resolves the production detector for a database (the active
  model, or the hybrid). `analyze_exposure!` is the only code path that
  detects peaks (the CLI's `analyze` and the ingest pass both go through it),
  so it is the one caller; code that only stamps versions, such as the
  migration that backfills `analysis_version`, must not stamp a model version.
- A promotion starts a **model pass** over unreviewed exposures whose
  `peak_model_version` differs, one at a time, newest experiments first, after
  any previous pass has finished. It calls `analyze_exposure!` with a
  `force_detect` flag, because the trace hash and analysis version are
  unchanged and would otherwise skip detection. Inside the write lock and
  transaction it re-checks `reviewed_at` and leaves the exposure alone if a
  curation landed meanwhile.
- The pass writes no `user_actions` rows but still sends each exposure's
  `analyze_run` frame. That frame is a cache update, not a notice: open
  clients replace stale auto peak ids, so a later exclude does not hit a
  deleted row. No assignment-dropped toast can follow, since reviewed
  exposures are never in the pass.
- After a core version bump, boot re-analysis of unreviewed exposures was
  designed to wait for the boot retrain (and its gate), so they go to the new
  model, or to the hybrid, once instead of hybrid then model. Step 4 does not
  wait (see "Model-written peaks (step 4)"). Reviewed exposures are not
  re-analysed on a bump.
- The model's output uses the pipeline's peak shape: indices on the trace
  grid, sharpness from `sharpness_full` with the hybrid's fallback (the median
  positive value for sharpness ≤ 0), prominence from persistence or 0 when the
  point is not a local maximum, and `found_by` passed to
  `diff_update_auto_peaks!`. If prediction fails the exposure gets the
  hybrid's peaks, not an analysis error.
- Prediction reuses the hybrid result (including `z`) that
  `analyze_exposure!` already computed instead of rerunning `findpeaks`.

## Model-written peaks (step 4)

What is built, and where it differs from "Analysis with a model" above.

**Resolver.** `production_detector(db)` returns the newest `active` row as a `ProductionModel` (id and
`TrainedModel`) or `nothing` for the hybrid. It reads the active id with one query per call and keeps the
parsed model in memory per database file, keyed by that id, so a promotion or retire (a new or missing active id) is
seen on the next call without a hook. Promotion and `invalidate_peak_model_cache!` empty it, which an edit to a row that keeps its id needs. The model is rebuilt
from the row: the recipe, the fitted learner and the stored feature version. A row whose feature version differs from
`feature_version` of its recipe's features, or that cannot be parsed, is not used (a warning is logged and the
hybrid is the detector). That is what keeps `bin/himalaya analyze` and the ingest pass from using a stale model after a core
upgrade, whatever order the server boots in.

**Analysis.** `analyze_exposure!` runs the hybrid, then, for an exposure whose `reviewed_at` is NULL, the model
over the same candidate pool (the hybrid result is passed in, so `findpeaks` runs once). The peaks use the pipeline shape: grid indices,
sharpness from the hybrid's `sharpness_full` with the median-positive fallback (the function the gate uses), prominence from `persistence`
or 0 off a local maximum. `found_by` is `both` when a hybrid peak is within `peak_match_tol` and `model` otherwise,
`auto_peaks.model_id` is the model and `exposures.peak_model_version` is set (NULL when the hybrid wrote the
peaks). Every detection now passes `found_by` explicitly, so a re-detection by the hybrid rewrites a row's provenance
to `hybrid`. Falls back to the hybrid, with a warning and never an analysis error, when prediction throws and when the trace
is shorter than the feature code's minimum (`MIN_POINTS`, where the candidate pool is empty). The decision is
checked again inside the write transaction: an exposure that was pinned meanwhile, or whose active model changed, gets the hybrid's peaks
(or, in the pass, is skipped). `peak_model_version` is not on the wire (the exposure list strips it like `reviewed_at`).

**Snapshots.** `record_review!` writes `shown_model_id` from the exposure's `peak_model_version` at event time.

**Promotion.** The server records a `passed` run and promotes it in the same `_DB_WRITE_LOCK` section
(`promote_peak_model!`, which parses the new model before taking the lock and seeds the resolver cache with it): when nothing usable is active, or when the run's `counts.curations` (add and exclude
rows, counted by the cold-start check at the run's own watermark) is at least 25% above the active run's
(`PEAK_MODEL_GROWTH`; an active run with no count is replaced by any run with one). The previous active row
becomes `retired` and loses its fitted JSON. A run whose model fails to load is left `passed`. Newest `passed` and
`active` keep their fitted JSON, as before. A passing run that is not promoted stays `passed` and is replaced
by the next `passed` run, so promotion compares each run only with the model in use.

**Model pass.** `start_peak_model_pass!` (after a promotion, and at boot when a model is active) starts one background task
(`PEAK_MODEL_PASS`: lock, `running`, `pending`, `stopped`). A start during a pass sets `pending`, and the pass
sweeps again when it ends, so a promotion during a sweep is never missed. A sweep lists exposures with `reviewed_at IS NULL`, a trace hash,
the current `ANALYSIS_VERSION`, no analysis error and `peak_model_version` different from the active id, newest experiment first,
and handles them one at a time on its own connection, yielding between exposures. Each goes through
`analyze_exposure!(...; force_detect = true, only_unreviewed = true, log_event = false)`. `force_detect` re-detects despite
matching hashes. `only_unreviewed` re-checks `reviewed_at` and the active id inside the transaction and skips the
exposure otherwise (so a curation that lands during the prediction is never overwritten, and a prediction failure leaves the
exposure for the next sweep rather than writing hybrid peaks). `log_event = false` writes no `user_actions` row. The pass sends no ingest frames. It does send the
`analyze_run` frame, with `id` 0 (no event row), so open clients replace stale peak ids. A sweep repeats (at most 5 rounds) until no stale exposure remains, which also catches an exposure an ingest wrote with the hybrid while a promotion landed. An exposure the model could not write (missing trace, a trace too short, a prediction that always throws) is remembered in memory per model id and not retried until the next promotion or restart; each run logs one summary warning for them instead of one per exposure. An error on one exposure does not stop the sweep, and the pass does not write `analysis_error`. An exposure whose peak set and index inputs come out unchanged sends no `analyze_run` frame, so a second promotion that changes nothing is silent. `stop_peak_model_pass!` sets `stopped` (cleared by
`start_background!` when it resumes work) and waits for the task. The frame needed a client change: `applyRemoteToCache` invalidated the
peaks cache of an `analyze_run` frame only when the trace hash changed, and now also when `findpeaks_skipped` is false (detection ran).

**Boot and retrain.** `start_background!` calls `boot_peak_model!` when it resumes work: it retires `active` rows whose
feature version is out of date (`retire_stale_peak_models!`) and starts a pass if a model is active. A retrain is due, with no new
label events and once per process, when the newest non-failed `peak_models` row has a feature version other than the current one
and labels exist (`peak_training_stale`).

Deviation: the design waited for the boot retrain before re-analysing unreviewed exposures after a core bump, so they go
to the new model once. This build does not wait. Boot analysis of stale exposures runs straight away with the
hybrid (the stale model is retired or unusable), and the retrain, which takes minutes in a child process, promotes a model that the pass then applies. Unreviewed
exposures are therefore analysed twice after a core bump, but ingest and the first view never wait on training, and a failed
or rejected retrain leaves the hybrid in place without any extra state.

**Cold start.** The precompile workload builds a three-tree model, promotes it and runs a pass over a synthetic trace, and fails
the build if the model did not write that exposure's peaks. Model JSON is read through `JSON.parse(::String)`: with an abstract
`AbstractString` the precompiled analysis code was thrown away at load (see "Cold start" in `packages/HimalayaUI/src/AGENTS.md`).

## No model UI

The model has no status view, controls or notices in the app. Each run's
outcome (promoted or not, and why, or the cold-start reason) is recorded in
its `peak_models` row and logged by the server, for whoever runs the
deployment.

## Packaging

Done in step 1. The features, labels, model and cross-validation live in
`packages/HimalayaPeakModel`, which depends only on Himalaya, DecisionTree and stdlib packages.
The database loader is `HimalayaUI.labeled_exposures(db)`, taking an open
handle. The command-line tool `tools/peak-classifier/` sits on top of both and
keeps its read-only open. HimalayaUI depends on HimalayaPeakModel, so
DecisionTree is in the sysimage; prediction goes into the precompile workload
in step 4.

## Privacy

Nothing leaves the machine: the model is trained from the local database and
trace files and stored in that database. The forest, metadata and feature
cache are lab data like the rest of the database and are never exported.

## Delivery

1. Packaging: move features, model and cross-validation into the new package;
   loader takes a database handle; the tool keeps its CLI, and excludes that
   no longer match an auto peak now count as negatives.
2. Provenance and pinning: `found_by`, `model_id`, `reviewed_at`, review
   snapshots and the snapshot-based labels. No peak changes.
3. Training service: `peak_models` table (recipe, recipe id and fitted model
   as JSON), feature cache, child-process
   training, watermark trigger and the gate. Runs are recorded but nothing is
   promoted, so peaks do not change. Done; see "Training service (step 3)".
4. Model-written auto peaks: `peak_model_version`, prediction from the stored JSON model,
   the model pass with `force_detect`, boot ordering after a core bump, and
   the one detector resolver for server and CLI. Done; see "Model-written peaks (step 4)".

## Evidence

Local measurements on a curated real database (leave-one-experiment-out, out
of fold), reported here in aggregate terms only.

- **Phase calls.** With classifier peaks scored through the same
  `auto_group` rule and sharpness as the hybrid's, top-1 phase calls on
  confirmed exposures tie the hybrid at the default threshold (with roughly
  equal paired wins and losses) and are worse at a stricter
  threshold. The gate therefore requires calls not to fall rather than to
  improve.
- **Missed peaks.** The classifier recovers markedly more of the peaks users
  added by hand than the hybrid does, with a few more peaks per exposure.
  Those adds were selected by what the older detector missed, so this measures
  recovery of that detector's misses. That is the expected gain: less manual
  peak work.
- **Detector dependence.** After the database is re-analysed with the hybrid,
  the hybrid scores perfectly on kept and excluded-detected by construction,
  and a large share of excludes no longer match any auto peak. Hence the
  hybrid-only comparison and the exclude labelling above.
- **Label era.** The labels behind these measurements were made while an older
  detector was on screen, neither the hybrid nor the classifier, so the
  comparison between those two is not biased by which one users saw.
