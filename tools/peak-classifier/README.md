# Peak classifier (local training tool)

Learns which candidate peaks are real from the curations in your own
HimalayaUI database: auto peaks you kept or excluded and peaks you added by
hand. It is not wired into the app. It is a separate Julia project, so the
core package and HimalayaUI gain no dependencies.

## Privacy

Everything stays on your machine. The tool opens the database with the SQLite
URI `mode=ro` (plus `PRAGMA query_only`), so it runs no migration, creates no
file for a mistyped path, and does not checkpoint a `-wal` file when it closes.
It reads your trace files and writes only the model path you give it. A
database in WAL mode (one the server has open) needs its `-shm` and `-wal`
files to exist or its directory to be writable, so a read-only open may leave
empty `-shm` and `-wal` side files next to it; those are the only files it may
add. The safest input is a copy of the database.

Trained models, feature files and evaluation outputs are lab data: the
repository ignores `tools/peak-classifier/models/`, `features/` and
`evaluations/`. Keep them there and never commit one.

## Run

From `tools/peak-classifier/`, once: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.

```
julia --threads=auto --project=. scripts/peak_classifier.jl cv    DB
julia --threads=auto --project=. scripts/peak_classifier.jl train DB --out models/lab.json
julia --project=. scripts/peak_classifier.jl predict models/lab.json TRACE_tot.dat
julia --threads=auto --project=. scripts/peak_classifier.jl evaluate DB --detector hybrid,persistence,recipe --json evaluations/run.json
```

Options for `cv` and `train`: `--learner rf|logistic` (default `rf`), `--threshold` (default 0.5), `--nms` (8 grid
points), `--trees` (200, forest only), `--group-by experiment|sample` (cross-validation
fold key), `--no-form-factor`, `--no-detector-features` (drop the
`in_persistence` and `in_hybrid` features; see "Reading the metrics").
`predict` prints one peak q per line. From Julia:
`predict_peaks(load_model(path), q, I, σ)`.

A model file is JSON: `{"file_version": 3, "model": ...}`, where `model` is
`HimalayaPeakModel.to_dict` of the trained model (recipe, recipe id, feature
names, threshold, NMS window and the fitted model, each tree as nested split and
leaf nodes). Loading reads data only and runs no code from the file, and it does
not depend on the Julia version (it is rebuilt into DecisionTree.jl 0.12's
tree structs). A file from an older version of the tool (a `Serialization`
dump, or file version below 3), or one trained with a different feature
version (for example after a core Himalaya upgrade), is refused with "retrain it". Any learner kind registered in `HimalayaPeakModel`
(`register_learner!`) can be saved and loaded; `--learner` selects among the
built-in `rf` (random forest) and `logistic`. Pass
the trace's `σ` to `predict_peaks` when you have it; without it `σ` is
estimated, and predictions can differ from what was validated.

Tests: `julia --project=. -e 'using Pkg; Pkg.test()'` (synthetic data only; CI job
"Peak classifier tool").

## Evaluate

`evaluate DB` scores peak detectors on one database, with the same scoring
code the server's training gate uses (`HimalayaUI.peak_eval`). The `hybrid` row
agrees with a recorded run's `comparison.hybrid`. A `recipe` row agrees with its
`comparison.model` only at the gate's chosen threshold (pass it as `--threshold`)
and with the same learner options; the default here is 0.5.

Options: `--detector` (comma list of `hybrid`, `persistence`, `recipe`, `oracle`; default
`hybrid`), `--threshold X|sweep` (recipe decision threshold; `sweep` scores
every threshold of the gate's grid, 0.3 to 0.9 and omits the recipe id; ignored with a warning when `recipe` is not a detector), `--learner`, `--nms`, `--trees`,
`--no-detector-features` (recipe, as for `cv`), `--phase-auto user_events|all`,
`--phase-targets reachable|all` (default `reachable`),
`--group-by experiment|sample` (fold key; by default training's choice,
experiments when there are at least 3 and otherwise samples; with fewer than
5 samples training refuses, and `evaluate` still runs, by sample), `--json PATH`.

Every loaded exposure of the database is scored, labelled or not. Output, per
detector, as numerator/denominator and rate:

- **phase top-1**: of the labelled exposures that have a reachable confirmed call, how many
  have the first pick of `auto_group` (computed from the detector's peaks, minus
  peaks the user excluded) matching a confirmed call in phase and basis (2%).
  A second table gives this per phase, using each exposure's first confirmed call.
  A call is reachable when `indexpeaks` gives a candidate of the called phase
  with a basis within 2% of the call, over the oracle's peaks or over the peaks
  of any detector in the run. So the targets depend on the detector set (the
  oracle is always included) and are the same for every detector in a run.
  Calls no peak set produces (for example custom indices anchored on fewer peaks
  than indexing needs) only dilute the metric and are left out; an exposure with
  no reachable call is not scored, and one with any keeps all its calls. Calls
  whose phase exists but is out-ranked stay in. `--phase-targets all` scores
  every call. The output gives reachable and all counts, overall and per phase.
- **add recall**: your added peaks the detector finds. Adds in the top 7% of
  the q range are not scored (the hybrid trims it); their count is printed.
- **excludes detected**: excluded peaks the detector still outputs (lower is better).
- **oracle** (a reference, not a detector): the user's own final peaks for each
  labelled exposure (kept auto peaks and adds, snapped to the nearest grid point,
  excludes removed), through the same scoring. Its add recall is close to 1, not
  exactly 1: adds within match tolerance of an exclude are dropped with it, and
  adds outside the trace's q range cannot be placed. Its phase top-1 is not an
  upper bound: a detector's different peak set can make `auto_group`'s first pick
  match a call the user's peaks miss. Its peaks/exposure counts only labelled
  exposures' peaks, so ignore that row.
- **peaks/exposure** over all loaded exposures, and **peaks/form-factor** over
  form-factor exposures.

With more than one detector, each is compared with the first on the same
exposures: per phase, wins (this detector right, the first wrong) and losses.
Wins minus losses is the difference in correct calls.

`recipe` is out of fold: a model trained on all other groups scores each held-out
group (leave one experiment out, or one sample with fewer than 3 experiments).
`hybrid` and `persistence` (`findpeaks` defaults) have no fitting, so the
whole database is held out by construction.

Biases to keep in mind:

- Labels and confirmed calls were made while some detector's peaks were on
  screen, so that detector is favoured. Only exposures last reviewed without a
  model count, and `--phase-auto all` also counts auto calls nobody edited, which
  favours whatever detector produced them.
- A confirmed call is a proxy for the true phase: it is the user's latest
  accepted call, and a wrong or loosely accepted call counts as truth.
- A threshold sweep picks the best row on the data it scores. Treat the best
  threshold's numbers as optimistic and compare detectors at a threshold fixed
  beforehand.
- With few experiments the folds are large and per-phase counts are small.
- The phase metric can stay below 100% even with the user's own peaks: it
  needs the first `auto_group` pick to match a confirmed call within 2% of its
  basis, and a reachable call can still be out-ranked (coexisting phases). Read
  each detector against the oracle row of the same phase. With `--phase-targets
  all`, unreachable calls also count against every detector.

The JSON holds the same numbers (`counts` including `phase_calls`, `phase_targets`, `detectors[].metrics`,
`detectors[].per_phase`, `paired`), with no paths. The output is derived from
lab data: keep it under `evaluations/` and do not commit it.

## How it works

1. **Candidates.** The union of the persistence detector, the hybrid detector
   and local maxima of the hybrid's matched-filter z above 3 (at least 3
   points apart, top 7% of the q range excluded as in the hybrid). This pool
   is deliberately loose: it sets the recall ceiling.
2. **Features** (`FEATURE_NAMES`). Matched-filter z per kernel width and the
   noise-normalised raw z, height above the SNIP background, persistence
   prominence, curvature, background slope, position, local density of
   candidates, and whether the persistence and hybrid detectors found it. All
   are noise-normalised or relative, so intensity units do not matter.
3. **Labels** (`label_candidates`). A labelled exposure is one where you
   edited its own peaks or indexing (a peak add, remove, exclude or un-exclude,
   or an assignment edit); phase calls on a sample and accept or reject do not
   count. Its
   kept peaks come from the latest direct review snapshot, the auto peaks you saw, minus
   those an exclude or add now hides. Rules apply in order, first match wins:
   - positive: the nearest candidate to each add, and to each kept peak the
     hybrid detector found (`found_by` hybrid or both);
   - dropped: other candidates inside the match tolerance of a positive
     (duplicates; non-maximum suppression handles them at prediction time),
     and candidates inside the tolerance of a kept peak only a model found;
   - strong negative: a candidate within the match tolerance of any exclude
     you recorded (even one no stored auto peak matches any more), and every
     remaining candidate on a form-factor exposure;
   - weak negative: any other candidate on a labelled exposure. Weak negatives
     are subsampled at 50% in training.
   Other exposures give no labels. Tolerance and the matching of curations to
   auto peaks come from `peak_match_tol` (defined in `HimalayaPeakModel`, used
   by HimalayaUI) and `HimalayaUI.curation_flags`, the app's own code. The
   database must have been opened once by a HimalayaUI that has review
   snapshots (`peak_reviews`); the tool does not migrate.
4. **Model.** A recipe (`HimalayaPeakModel.Recipe`): feature set, learner,
   label options and decision rule. The default is the random forest
   (DecisionTree.jl); `--learner logistic` swaps in a logistic regression.
   Prediction keeps candidates with
   probability at or above the threshold, then thins them by non-maximum
   suppression.

## Reading the metrics

`cv` holds out one group at a time (an experiment by default), trains on the
rest and scores the held-out group, so every number is out of fold. Rows are
thresholds, the hybrid detector, and the unfiltered candidate pool. Columns:

- **add recall**: fraction of your added peaks the output contains. The
  all-candidates row is the ceiling.
- **kept recall**: fraction of kept auto peaks still present.
- **excl detected**: fraction of excluded auto peaks still present. Lower is
  better.

Kept recall and excl detected depend on which detector wrote `auto_peaks`. On a
database analysed with the hybrid detector, the hybrid row has kept recall and
excl detected near 1.0 by construction, and the `in_hybrid` and
`in_persistence` features line up with the kept and excluded labels. Treat add
recall and extra/labelled as the independent columns. `--no-detector-features`
trains without those two features, to check how much the model leans on them.
Other columns:

- **extra/labelled**: predicted peaks per labelled exposure that match no add and
  no kept peak. Real peaks you never labelled count as extras, so this
  overstates false positives.
- **peaks/FF**: predicted peaks per form-factor exposure.
- **candidate AUC**: ranking quality of the probabilities, positives against
  all negatives.

Pick the threshold from this table. The default of 0.5 came from a sweep on
one other dataset; choosing from the same table you evaluate on is slightly
optimistic. Every experiment is a different group, so with few experiments the folds are large
and noisy. Compare against the hybrid row, not against 1.0.

## Label pitfalls

- Users sometimes exclude real broad sample peaks, so some negatives are real
  peaks. The model then learns that inconsistency.
- Unlabelled candidates on a curated exposure may be real peaks you did not
  bother to add (weak negatives).
- Kept auto peaks are implicit accepts. They are only trusted on exposures
  with at least one curation, since an exposure you never opened gives no signal.
- Stored auto peaks come from whichever detector version ran when the data
  was analysed.
- Instrument artefacts (a detector-gap step, a beamstop edge) produce
  confident false candidates. Marking such exposures form factor, or
  excluding the artefact peaks, teaches the model to reject them.
