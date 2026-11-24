using DBInterface, SQLite

# physics: I_2 can collapse to a few percent of normal while I_0 stays normal, for
# physics: example when the beam misses the beamstop diode. SASTool still divides that
# physics: image's trace by its I_2, which inflates the trace many times over, and the
# physics: beam-off test does not catch it. An image whose I_2 per I_0 per counting
# physics: second is below a quarter of its setup's median counts as collapsed.
const I2_COLLAPSE_FRACTION = 0.25

_positive_reading(v) = !ismissing(v) && v !== nothing && isfinite(Float64(v)) && Float64(v) > 0

transmission_per_flux(transmitted, incident, counting_time) =
    _positive_reading(transmitted) && _positive_reading(incident) && _positive_reading(counting_time) ?
        Float64(transmitted) / (Float64(incident) * Float64(counting_time)) : nothing

function setup_transmission_median(db::SQLite.DB, setup_id)::Union{Float64,Nothing}
    (ismissing(setup_id) || setup_id === nothing) && return nothing
    rows = sqlrows(DBInterface.execute(db,
        "SELECT transmitted_counts, incident_monitor, counting_time FROM exposures
          WHERE setup_id = ? AND kind = 'file'
            AND transmitted_counts > 0 AND incident_monitor > 0 AND counting_time > 0",
        [Int(setup_id)]))
    values = sort!(Float64[t for t in (transmission_per_flux(r.transmitted_counts, r.incident_monitor, r.counting_time) for r in rows) if t !== nothing])
    isempty(values) && return nothing
    n = length(values)
    isodd(n) ? values[(n + 1) ÷ 2] : (values[n ÷ 2] + values[n ÷ 2 + 1]) / 2
end

function is_i2_collapsed(transmitted, incident, counting_time, median)::Bool
    median === nothing && return false
    _positive_reading(incident) && _positive_reading(counting_time) || return false
    (ismissing(transmitted) || transmitted === nothing) && return false
    t = _positive_reading(transmitted) ? Float64(transmitted) / (Float64(incident) * Float64(counting_time)) : 0.0
    t < I2_COLLAPSE_FRACTION * median
end
