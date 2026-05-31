# Gauss–Bonnet lattice-parameter ratios among the three bicontinuous cubic
# phases. Reference scale (Pn3m = 1.0): a_Pn3m : a_Im3m : a_Ia3d.
#
# The 1.000:1.279:1.576 triple is the canonical lattice ratio for the
# D (Pn3m) / P (Im3m) / G (Ia3d) infinite periodic minimal surfaces related by
# the Bonnet (isometric) transformation — standard in the bicontinuous-cubic-
# phase literature (e.g. Hyde et al., "The Language of Shape", 1997). It follows
# from the Gauss–Bonnet relation between each surface's area and Euler
# characteristic. NOTE: this is NOT the ratio of the NGC A₀ surface-area
# constants in index.jl (1.919/2.345/3.091) — those differ (2.345/1.919 = 1.222,
# not 1.279). `bonnet_lattice` returns `nothing` for any non-bicontinuous-cubic
# pairing — the Bonnet relation is cubic↔cubic only.
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

"""
    bonnet_consistent(from, a_from, to, a_to_observed; reltol=0.02) -> Bool

True when an observed coexisting-cubic lattice `a_to_observed` matches the
Bonnet-predicted lattice within relative tolerance `reltol`. False for phase
pairs with no Bonnet relation.
"""
function bonnet_consistent(from::Type{<:Phase}, a_from::Real,
                           to::Type{<:Phase}, a_to_observed::Real; reltol::Real=0.02)
    pred = bonnet_lattice(from, a_from, to)
    pred === nothing && return false
    abs(a_to_observed - pred) <= reltol * pred
end
