# Himalaya.jl
SAXS diffraction pattern indexing the easy way.

This repository contains:

- **`Himalaya`** — the core Julia library for peak-finding and phase indexing (documented below).
- **`HimalayaUI`** — a web application that wraps the core library with an interactive UI for batch-analyzing a beamtime experiment and refining phase assignments. See [Web application](#web-application).

## Installation
Once you have started Julia, enter [Pkg mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Pkg-mode) by typing `]` and run `add Himalaya` command. You will see something like this
```julia-repl
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.8.0-beta3 (2022-03-29)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

(@v1.8) pkg> add Himalaya
```
And that's it!

## Usage
The following code is an example of how to use Himalaya to index the integration
of a diffraction pattern.

```julia
using DelimitedFiles
using Himalaya

# `_tot.dat` files are space-separated `q  I(q)  σ(q)` (third column is the
# per-point intensity uncertainty from azimuthal integration).
A = readdlm("my-high-impact-sample_tot.dat", ' ', Float64, '\n')
qs, Is, σs = A[:, 1], A[:, 2], A[:, 3]

# Detect peaks using topological prominence + curvature, with adaptive
# (kneedle) thresholds on both. Shape-agnostic; no per-trace tuning.
pk = findpeaks(qs, Is, σs)        # NamedTuple: (indices, q, prominence, sharpness)
peak_qs    = pk.q
peak_proms = pk.prominence

# Index the detected peaks against known phases.
indices = indexpeaks(peak_qs, peak_proms)
```

You can also get `Index`s for a specific phase as follows
```julia
indexpeaks(Hexagonal, peak_qs, peak_proms)
# 6-element Vector{Index}:
#  Index(::Hexagonal, 0.06426, [0.06426 0.11144 0.12876 0.16771 0.19368 ⋅ ⋅])
#  Index(::Hexagonal, 0.07855, [0.07855 0.13655 0.15732 ⋅ ⋅ ⋅ ⋅])
#  Index(::Hexagonal, 0.09673, [0.09673 0.16771 0.19368 0.25645 0.29107 ⋅ ⋅])
#  ...
```
You can use `?Phase` to see the list of available phases (`?` enters [help mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Help-mode)).

## Web application

`HimalayaUI` is a React + Julia (Oxygen.jl) web app that runs against a beamtime experiment directory. Point it at the experiment folder and a sample manifest, and it auto-indexes every exposure and presents the results in a browser for review and refinement.

```bash
# 1. Initialize a SQLite DB from the manifest and experiment layout
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  init /path/to/experiment --manifest manifest.csv

# 2. Run batch analysis (peak-finding + indexing + auto-grouping, all samples)
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  analyze /path/to/experiment

# 3. Start the web UI (default port 8080; open http://localhost:8080)
julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' -- \
  serve /path/to/experiment --port 8080
```

The UI shows a trace viewer (log-log I(q) with σ ribbon), a Miller-index scatter with linear fit, and a phase-assignment panel with per-exposure auto-picked and alternative candidates. Peaks can be added or removed by clicking the trace; indices can be confirmed or excluded from the active group. All edits are attributed and audited via a username header.

See [`packages/HimalayaUI/`](packages/HimalayaUI/) for the sub-package, and [`docs/superpowers/specs/`](docs/superpowers/specs/) for the full design spec.
