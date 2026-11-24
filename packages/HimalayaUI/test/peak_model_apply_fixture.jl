using Test, HTTP, Random, SQLite, DBInterface, Tables
import JSON
using HimalayaUI
import HimalayaPeakModel as HPM
isdefined(@__MODULE__, :_pt_db) || include("peak_training_fixture.jl")

_pma_req() = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
_pma_one(db, sql, args = []) = only(Tables.rowtable(DBInterface.execute(db, sql, args)))

const _PMA_TRAINED = Ref{Any}(nothing)
function _pma_trained()
    _PMA_TRAINED[] === nothing || return _PMA_TRAINED[]
    mktempdir() do dir
        b = _pt_db(dir)
        DBInterface.close!(b.db)
        ro = HimalayaUI.open_readonly_db(b.path)
        r = HimalayaUI.run_peak_training(ro; recipe = HPM.default_recipe(n_trees = 15), io = devnull,
                                         peaks_slack = 100, ff_slack = 100)
        DBInterface.close!(ro)
        _PMA_TRAINED[] = JSON.parse(JSON.json(r))
    end
    _PMA_TRAINED[]
end

function _pma_result(; curations = nothing, watermark = nothing)
    r = Dict{String, Any}(String(k) => v for (k, v) in pairs(_pma_trained()))
    r["status"] = "passed"
    curations === nothing || (r["counts"] = merge(Dict{String, Any}(String(k) => v for (k, v) in pairs(r["counts"])),
                                                  Dict("curations" => curations)))
    watermark === nothing || (r["watermark"] = watermark)
    r
end

function _pma_add_model!(db; status = "passed", kw...)
    id = HimalayaUI.insert_peak_model!(db, _pma_result(; kw...))
    status == "passed" || DBInterface.execute(db, "UPDATE peak_models SET status = ? WHERE id = ?", [status, id])
    id
end

_pma_status(db, id) = _pma_one(db, "SELECT status, fitted_json IS NOT NULL AS has_fitted FROM peak_models WHERE id = ?", [id])

function _pma_db(dir; n = 3, points = 700, experiments = 1)
    rng = Random.MersenneTwister(5)
    db = open_prepared_clone(dir)
    xids = Int[]; sids = Int[]; eids = Int[]; adirs = Dict{Int, String}(); traces = Dict{Int, Any}()
    for x in 1:experiments
        adir = joinpath(dir, "analysis$x"); mkpath(adir)
        xid = HimalayaUI.init_experiment!(db; path = joinpath(dir, "exp$x"), data_dir = joinpath(dir, "data$x"), analysis_dir = adir)
        sid = HimalayaUI.create_sample!(db; experiment_id = xid, name = "S$x")
        push!(xids, xid); push!(sids, sid)
        for k in 1:n
            name = "x$(x)k$(k)_tot"
            q, I, σ = _pt_trace(rng, 0.05 + 0.01k; n = points)
            _pt_write_dat(joinpath(adir, name * ".dat"), q, I, σ)
            eid = HimalayaUI.create_exposure!(db; experiment_id = xid, sample_id = sid, filename = name)
            push!(eids, eid); adirs[eid] = adir; traces[eid] = HimalayaUI.load_dat(joinpath(adir, name * ".dat"))
        end
    end
    (; db, adir = adirs[first(eids)], adirs, xid = first(xids), xids, sid = first(sids), sids, eids, traces)
end

_pma_peaks(db, e) = Tables.rowtable(DBInterface.execute(db,
    "SELECT q, found_by, model_id, sharpness, prominence, findpeaks_index FROM auto_peaks WHERE exposure_id = ? ORDER BY q", [e]))
_pma_version(db, e) = _pma_one(db, "SELECT peak_model_version AS v, analysis_version FROM exposures WHERE id = ?", [e])
_pma_review!(db, e) = DBInterface.execute(db, "UPDATE exposures SET reviewed_at = '2026-01-01T00:00:00.000Z' WHERE id = ?", [e])

