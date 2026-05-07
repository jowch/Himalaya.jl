# src/ — Core Himalaya Package

## OVERVIEW
Peak-finding and phase-indexing algorithms for 1D SAXS traces.

## WHERE TO LOOK
| Task | Location | Notes |
|------|----------|-------|
| Peak finding | `peakfinding.jl` | AND-gate of prominence + sharpness; adaptive kneedle thresholding |
| Persistence helper | `persistence.jl` | Topological persistence via Peaks.jl |
| Sharpness | `sharpness.jl` | Savitzky-Golay 2nd derivative or CWT |
| Thresholding | `threshold.jl` | Kneedle elbow finder |
| Phase types | `phase.jl` | Phase abstract hierarchy + phaseratios |
| Indexing | `index.jl` | Index struct, indexpeaks, score |
| Utilities | `util.jl` | Shared math helpers |

## CONVENTIONS
- TDD: failing test → minimal impl → verify → commit
- Regression floors, not hard-coded counts (`recall ≥ floor`, `spurious ≤ ceiling`)
- `include()` order must respect dependency graph (no forward refs)

## ANTI-PATTERNS
- `Fm3m` is missing from `indexpeaks` dispatch — known gap, don't fix opportunistically
- `score(index)` = `coverage × consistency`; R² is NOT part of score
