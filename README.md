# Himalaya.jl
SAXS diffraction pattern indexing the easy way.

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

load_integration(path) = readdlm(path, ' ', Float64, '\n')

# `integration` contains the values in a tot_file
integration = load_integration("my-high-impact-sample_tot.dat")
qs, logIs = integration[:, 1], log10.(integration[:, 2])

# indices of peaks in the integration array
peak_locations, peak_proms = findpeaks(logIs)
peak_qs = qs[peak_locations]

# compute phases matching identified peaks
indices = indexpeaks(peak_qs, peak_proms)

# Example `indices`
# 3-element Vector{Index}:
#  Index(::Pn3m, 0.05344, [0.05344 0.06556 0.07508 0.09239 0.10668 0.11317 ⋅ 0.12486]
#  Index(::Im3m, 0.05344, [0.05344 0.07508 0.09239 0.10668 ⋅ 0.13049 ⋅ ⋅ ⋅]
#  Index(::Pn3m, 0.07508, [0.07508 0.09239 0.10668 0.13049 ⋅ ⋅ ⋅ ⋅]

index = first(indices)

score(index) # => 6.4996...
d, R² = fit(index) # => d = 117.8585...Å; R² = 0.9999...
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
