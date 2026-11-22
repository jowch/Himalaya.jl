# Peak classifier (local training tool)

Learns which candidate peaks are real from the curations in your own
HimalayaUI database: auto peaks you kept or excluded and peaks you added by
hand. It is not wired into the app. It is a separate Julia project, so the
core package and HimalayaUI gain no dependencies.

## Privacy

Everything stays on your machine. The tool reads the database file read-only
(no migration, `PRAGMA query_only`) and your trace files, and writes only the
model path you give it. Trained models and feature files are lab data: the
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
fold key), `--no-form-factor`. `predict` prints one peak q per line. From Julia:
`predict_peaks(load_model(path), q, I, σ)`.

A model file is a Julia `Serialization` dump: it loads only on the same Julia
and DecisionTree.jl versions. Retrain after upgrading.

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
3. **Labels** (`label_candidates`). On a curated exposure (one with at least
   one add or exclude):
   - positive: the nearest candidate to each add and to each auto peak you
     did not exclude;
   - dropped: other candidates inside the match tolerance of a positive
     (duplicates; non-maximum suppression handles them at prediction time);
   - strong negative: a candidate at an excluded auto peak, and every
     candidate on a form-factor exposure;
   - weak negative: any other candidate on a curated exposure (you looked and
     did not add it). Weak negatives are subsampled at 50% in training.
   Exposures without curations give no labels. Tolerance and the matching of
   curations to auto peaks come from `HimalayaUI.peak_match_tol` and
   `curation_flags`, the app's own code.
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
- **extra/curated**: predicted peaks per curated exposure that match no add and
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
