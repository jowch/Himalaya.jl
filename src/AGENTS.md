# src/ — Core Himalaya Package

Peak-finding and phase-indexing algorithms for 1D SAXS traces.

## Where to look

| Task | Location | Notes |
|------|----------|-------|
| Peak finding | `peakfinding.jl` | AND-gate of prominence + sharpness; adaptive kneedle thresholding |
| Persistence helper | `persistence.jl` | Topological persistence via Peaks.jl |
| Sharpness | `sharpness.jl` | Savitzky-Golay 2nd derivative or CWT |
| Thresholding | `threshold.jl` | Kneedle elbow finder |
| Phase types | `phase.jl` | Phase abstract hierarchy + `phaseratios` |
| Indexing | `index.jl` | `Index` struct, `indexpeaks`, `score` |
| Utilities | `util.jl` | Shared math helpers |

Background reading: `docs/peak-finding.md` (load-bearing), `docs/scoring.md`.

## Conventions

- TDD: failing test → minimal impl → verify → commit
- Regression floors, not hard-coded counts (`recall ≥ floor`, `spurious ≤ ceiling`)
- `include()` order must respect dependency graph (no forward refs)

## Index scoring

`score(index)` returns a value in `[0, 1]` — product of:

- **coverage**: harmonic-weighted fraction of expected peaks found, `1/rank` weight per position.
- **consistency**: `1/(1+CV)` of peak sharpnesses. Guard `cv` against zero mean before dividing (all-zero sharpness is valid and should score as consistent).

`totalprom` and the `prom` field on `Index` no longer exist — the struct now has `sharpness::SparseVector`.

`auto_group` and `remove_subsets` in `pipeline.jl` both depend on `score` ordering — correctness of auto-analysis flows from score quality.

R² is stored per index but is **NOT** part of `score`; it is a UI hard gate (threshold 0.98 in `PhasePanel`).

## Anti-patterns

- **`Fm3m` is missing from `indexpeaks` dispatch.** The all-phases loop in `index.jl` omits `Fm3m`. The phase is defined and `minpeaks`/`phaseratios` exist, but `indexpeaks` can never return an `Fm3m` index. Known pre-existing gap — don't fix opportunistically.
