# Per-database peak model

Status: design. Steps 1 (packaging) and 2 (provenance and pinning) are implemented: the features, labels,
model and cross-validation live in `packages/HimalayaPeakModel`, the database
loader is `HimalayaUI.labeled_exposures` (reading review snapshots), and `tools/peak-classifier/` (#486)
is the command-line tool on top of both. This document describes how HimalayaUI trains such a model for its own
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
remove or state change; a sample call (which pins every exposure of that
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

A core `ANALYSIS_VERSION` bump re-detects every exposure, pinned or not, as
today, with whatever detector is in production at that moment; assignments
can be dropped exactly as on any release today. Labels are unaffected, because
training reads review snapshots, not live auto peaks (below).

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

- Table `peak_models` (migration): id, created_at, status (`active`,
  `pending`, `rejected`, `retired`), label watermark, counts, threshold and
  non-maximum-suppression window,
  metrics JSON (the cross-validated tables and the production comparison),
  feature version, core version, the recipe (its `to_dict` as JSON), the
  recipe id, and the fitted model.
- The fitted model is `to_dict` of the learner's `Fitted` object, stored as JSON
  in a TEXT column, only for the active and a pending model; it is set to NULL
  on rejected and retired rows. For the random forest each tree is nested split
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
- After a core version bump, boot re-analysis of unreviewed exposures waits
  for the boot retrain (and its gate) to finish, so they go to the new model,
  or to the hybrid, once instead of hybrid then model. Reviewed exposures are
  re-analysed straight away with the hybrid.
- The model's output uses the pipeline's peak shape: indices on the trace
  grid, sharpness from `sharpness_full` with the hybrid's fallback (the median
  positive value for sharpness ≤ 0), prominence from persistence or 0 when the
  point is not a local maximum, and `found_by` passed to
  `diff_update_auto_peaks!`. If prediction fails the exposure gets the
  hybrid's peaks, not an analysis error.
- Prediction reuses the hybrid result (including `z`) that
  `analyze_exposure!` already computed instead of rerunning `findpeaks`.

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
   promoted, so peaks do not change.
4. Model-written auto peaks: `peak_model_version`, prediction from the stored JSON model,
   the model pass with `force_detect`, boot ordering after a core bump, and
   the one detector resolver for server and CLI.

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
