# Future Feature Ideas

Ideas that are intentionally out of scope for current development but worth
preserving for later planning.

## Extended lattice types

Monoclinic and tetragonal lattice indexing. These require fitting a lattice
parameter *vector* (2+ free parameters) rather than a single basis, which
changes how `indexpeaks` and `score` work internally. Will require extending
the `Phase` type hierarchy and the indexing engine. Design the extension points
when the need is concrete.

## Cross-experiment sample comparison

Query `sample_tags` and analysis results across multiple experiment databases
to compare samples with similar compositions (e.g., same lipid/peptide system
across beamtimes). The SQLite-per-experiment schema is migration-friendly:
a thin aggregation layer can open multiple DBs and JOIN across them via
SQLite's ATTACH mechanism.

## Auto-best-exposure selection

Heuristics for automatically selecting the best exposure for a given sample:
- Detect flares (anomalous low-q intensity spikes)
- Score exposures by signal quality (peak prominence, background level)
- Flag poor-signal exposures automatically

## Summary table page (web UI)

Full-screen tabular view of all samples in an experiment: confirmed phases,
lattice parameters, R², score. Sortable, filterable by tag. Export to CSV/JSON.

## Stacked / waterfall comparison page (web UI)

Multi-sample trace overlay with configurable I-offset between traces.
Publication-quality SVG export. Useful for visualizing phase transitions
across a sample series.

## Sub-pixel peak positions

Parabolic interpolation around detected peak maxima for more precise q values.
Currently peaks are returned at grid positions.

## Background subtraction in pipeline

Automated SNIP or asymmetric least-squares background subtraction as a
pre-processing step, particularly useful for traces with steeply falling
backgrounds. Currently handled manually (background-subtracted exposures are
ingested as derived exposures in the DB).
