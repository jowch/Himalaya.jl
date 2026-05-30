# Gauss–Bonnet lattice-parameter ratios among the three bicontinuous cubic
# phases. Reference scale (Pn3m = 1.0): a_Pn3m : a_Im3m : a_Ia3d.
#
# These track the NGC constants already hard-coded in index.jl (A₀ =
# 1.919/2.345/3.091), which derive from the same IPMS geometry. The
# Pn3m:Im3m:Ia3d = 1.000:1.279:1.576 triple is confirmed canonical
# (saxs-physics-reviewer, 2026-05-30). `bonnet_lattice` returns `nothing` for
# any non-bicontinuous-cubic pairing — the Bonnet relation is cubic↔cubic only.
const BONNET_SCALE = Dict{DataType,Float64}(
    Pn3m => 1.000,
    Im3m => 1.279,
    Ia3d => 1.576,
)

"""
    bonnet_lattice(from::Type{<:Phase}, a_from::Real, to::Type{<:Phase}) -> Union{Float64,Nothing}

Predict the lattice parameter a coexisting cubic of phase `to` should have,
given an assigned cubic of phase `from` with lattice `a_from`, via the
Gauss–Bonnet ratios. Returns `nothing` when either phase is not one of the
three bicontinuous cubics (no Bonnet relation defined).
"""
function bonnet_lattice(from::Type{<:Phase}, a_from::Real, to::Type{<:Phase})
    (haskey(BONNET_SCALE, from) && haskey(BONNET_SCALE, to)) || return nothing
    Float64(a_from) * (BONNET_SCALE[to] / BONNET_SCALE[from])
end
