# Himalaya.jl
SAXS diffraction pattern indexing the easy way.

## Usage
The following code is an example of how to use Himalaya to index the integration
of a diffraction pattern.

```julia
using Himalaya

# `integration` contains the values in a tot_file
qs, logIs = integration[:, 1], log10.(integration[:, 2])

# indices of peaks in the integration array
peak_locations = findpeaks(logIs)
peak_qs = qs[peak_locations]

# compute phases matching identified peaks
indices = indexpeaks(peak_qs)

# Example `indices`
# 3-element Vector{Index}:
#  Index(::Pn3m, 0.05344, [0.05344 0.06556 0.07508 0.09239 0.10668 0.11317 ⋅ 0.12486]
#  Index(::Im3m, 0.05344, [0.05344 0.07508 0.09239 0.10668 ⋅ 0.13049 ⋅ ⋅ ⋅]
#  Index(::Pn3m, 0.07508, [0.07508 0.09239 0.10668 0.13049 ⋅ ⋅ ⋅ ⋅]

index = first(indices)

score(index) # => 6.4996...
d, R² = fit(index) # => d = 117.8585...Å; R² = 0.9999...
```
