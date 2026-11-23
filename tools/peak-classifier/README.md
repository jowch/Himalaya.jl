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

Trained models and feature files are lab data: the
repository ignores `tools/peak-classifier/models/` and `features/`. Keep
models there and never commit one.

## Run

From `tools/peak-classifier/`, once: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.

```
julia --threads=auto --project=. scripts/peak_classifier.jl cv    DB
julia --threads=auto --project=. scripts/peak_classifier.jl train DB --out models/lab.jls
julia --project=. scripts/peak_classifier.jl predict models/lab.jls TRACE_tot.dat
```

Options for `cv` and `train`: `--threshold` (default 0.5), `--nms` (8 grid
points), `--trees` (200), `--group-by experiment|sample` (cross-validation
fold key), `--no-form-factor`, `--no-detector-features` (drop the
`in_persistence` and `in_hybrid` features; see "Reading the metrics").
`predict` prints one peak q per line. From Julia:
`predict_peaks(load_model(path), q, I, σ)`.

A model file is a Julia `Serialization` dump: it loads only on the same Julia
and DecisionTree.jl versions, and only with the version of this tool that wrote it (model files from before the HimalayaPeakModel package cannot be read). Retrain after upgrading. Loading it runs
arbitrary code from the file, so load only models you trained yourself. Pass
the trace's `σ` to `predict_peaks` when you have it; without it `σ` is
estimated, and predictions can differ from what was validated.

Tests: `julia --project=. -e 'using Pkg; Pkg.test()'` (synthetic data only; CI job
"Peak classifier tool").

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
4. **Model.** Random forest (DecisionTree.jl). Prediction keeps candidates with
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
