const PEAK_MATCH_POINTS = 6

"""
Tolerance for attaching a curation to a peak: 6 grid points when the grid step
is known, else `max(1e-6, 0.001·|q|)`. HimalayaUI uses this same function for
curation matching, so labels and the app agree.
"""
function peak_match_tol(q::Real, q_step)
    base = max(1e-6, 0.001 * abs(Float64(q)))
    (q_step === nothing || q_step === missing) && return base
    max(base, PEAK_MATCH_POINTS * Float64(q_step))
end

"""Median spacing of a q grid, or `nothing` if it has no positive finite spacing."""
function grid_step(q::AbstractVector{<:Real})
    length(q) >= 2 || return nothing
    d = sort!(Float64.(diff(q)))
    m = isodd(length(d)) ? d[(end + 1) ÷ 2] : (d[end ÷ 2] + d[end ÷ 2 + 1]) / 2
    (isfinite(m) && m > 0) ? m : nothing
end
