const PEAK_MODEL_PREDICT = Ref{Function}(HimalayaPeakModel.predict_peaks)

function model_peaks_result(pm::ProductionModel, q::Vector{Float64}, I::Vector{Float64}, σ::Vector{Float64},
                            hybrid::NamedTuple, q_step::Union{Float64, Nothing})
    length(I) < HimalayaPeakModel.MIN_POINTS && return nothing
    qs = PEAK_MODEL_PREDICT[](pm.model, q, I, σ; hybrid)
    idx = Int[findfirst(==(v), q)::Int for v in qs]
    persistence = Himalaya.persistence(I)
    prominence_at = Dict(zip(persistence.indices, persistence.prominence))
    found_by = String[any(h -> abs(q[i] - h) <= peak_match_tol(q[i], q_step), hybrid.q) ? "both" : "model"
                      for i in idx]
    (result = (q = q[idx], indices = idx, prominence = Float64[get(prominence_at, i, 0.0) for i in idx],
               sharpness = _sharpness_at(hybrid.sharpness_full, idx), sharpness_full = hybrid.sharpness_full),
     found_by = found_by)
end

peak_reviewed(db::SQLite.DB, exposure_id::Integer) =
    !isempty(sqlrows(DBInterface.execute(db,
        "SELECT 1 AS x FROM exposures WHERE id = ? AND reviewed_at IS NOT NULL", [Int(exposure_id)])))

active_peak_model_id(db::SQLite.DB) = begin
    rows = sqlrows(DBInterface.execute(db,
        "SELECT id FROM peak_models WHERE status = 'active' ORDER BY id DESC LIMIT 1"))
    isempty(rows) ? nothing : Int(rows[1].id)
end

function predict_model_peaks(db::SQLite.DB, exposure_id::Integer, q, I, σ, hybrid, q_step)
    pm = production_detector(db)
    pm === nothing && return (nothing, nothing)
    peak_reviewed(db, exposure_id) && return (pm, nothing)
    try
        (pm, model_peaks_result(pm, q, I, σ, hybrid, q_step))
    catch err
        @warn "peak model $(pm.id) failed on exposure $exposure_id; using the hybrid detector's peaks" exception = err
        (pm, nothing)
    end
end
