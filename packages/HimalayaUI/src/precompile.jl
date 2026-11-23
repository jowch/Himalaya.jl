using PrecompileTools: @setup_workload, @compile_workload
using ImageCore: Gray
using ImageCore.FixedPointNumbers: Q0f31

const _WORKLOAD_PHASES = (Himalaya.Lamellar, Himalaya.Hexagonal, Himalaya.Square, Himalaya.Pn3m,
                          Himalaya.Im3m, Himalaya.Ia3d, Himalaya.Fm3m, Himalaya.Fd3m)

const _WORKLOAD_MIXED_BASES = (0.071, 0.113, 0.131, 0.061, 0.047, 0.039, 0.052, 0.033)

_workload_centers() = vcat(0.05 .* Himalaya.phaseratios(Himalaya.Pn3m)[1:6],
                           0.09 .* Himalaya.phaseratios(Himalaya.Hexagonal)[1:4])

function _workload_mixed_centers()
    centers = Float64[]
    for (P, b) in zip(_WORKLOAD_PHASES, _WORKLOAD_MIXED_BASES)
        ratios = Himalaya.phaseratios(P)
        append!(centers, b .* ratios[1:min(end, 6)])
    end
    filter!(<(0.45), sort!(unique!(round.(centers; digits = 4))))
end

function _workload_trace(path::AbstractString, centers = _workload_centers())
    q = collect(range(0.01, 0.5; length = 920))
    I = @. 2.0e3 * q^-1.5 + 50.0
    for c in centers
        @. I += 3.0e3 * exp(-((q - c) / 0.002)^2 / 2)
    end
    @. I += 0.3 * sqrt(I) * sin(7919.0 * q)
    open(path, "w") do io
        for k in eachindex(q)
            println(io, q[k], " ", I[k], " ", sqrt(I[k]))
        end
    end
end

function _workload_frame(path::AbstractString)
    n = 160
    counts = Int32[round(Int32, 2.0e3 * exp(-hypot(i - n / 2, j - n / 2) / 8)) for i in 1:n, j in 1:n]
    counts[1:4, :] .= -1
    FileIO.save(path, Gray.(reinterpret.(Q0f31, counts)))
end

function _workload_prp(path::AbstractString)
    write(path, """
        Image file name: frame.tif
        I_signal=1000.0 (corrected for dark counts)
        I_0(daqIon3.0)=1.5
        Exposure time=1
        Counting time=1.05
        Beam energy=9000.0 eV
        Pipe length=1700 mm
        Horizontal position=1.0 mm
        Detector=Pilatus 1M
        """)
end

function _workload_headers()
    ["Content-Type" => "application/json",
     "X-Username" => "precompile",
     "X-Client-Id" => "precompile",
     "X-Client-Op-Id" => string(rand(UInt64))]
end

function _workload_inproc(method::AbstractString, target::AbstractString, payload::Vector{UInt8})
    Oxygen.internalrequest(HTTP.Request(method, target, _workload_headers(), payload);
                           middleware = [db_connection_middleware], serialize = true, catch_errors = true)
end

function _workload_socket(port::Int)
    (method, target, payload) -> HTTP.request(method, "http://127.0.0.1:$port$target",
                                              _workload_headers(), payload; status_exception = false)
end

_recording(send, ::Nothing) = send
_recording(send, record::AbstractVector) =
    (method, target, payload) -> (push!(record, (method, target)); send(method, target, payload))

function _workload_call(send, method::AbstractString, target::AbstractString, body = nothing)
    payload = body === nothing ? UInt8[] : Vector{UInt8}(JSON3.write(body))
    res = send(method, target, payload)
    res.status < 300 || error("precompile workload: $method $target returned $(res.status)")
    json_body = occursin("json", HTTP.header(res, "Content-Type"))
    isempty(res.body) || !json_body ? nothing : JSON3.read(res.body)
end

function _workload_requests(send, experiment::Int, sample::Int, other::Int, exposure::Int, mixed::Int)
    _workload_call(send, "GET", "/api/health")
    _workload_call(send, "GET", "/api/users")
    _workload_call(send, "GET", "/api/experiments")
    _workload_call(send, "GET", "/api/experiments/$experiment")
    _workload_call(send, "GET", "/api/experiments/$experiment/bridges")
    _workload_call(send, "GET", "/api/experiments/$experiment/problems")
    _workload_call(send, "GET", "/api/samples")
    _workload_call(send, "GET", "/api/samples?experiment_id=$experiment")
    _workload_call(send, "GET", "/api/samples/$sample/exposures")
    _workload_call(send, "POST", "/api/samples/$sample/bridge", Dict(:to_sample_id => other))
    _workload_call(send, "DELETE", "/api/samples/$sample/bridge")
    _workload_call(send, "GET", "/api/experiments/$experiment/compares")
    _workload_call(send, "PATCH", "/api/samples/$sample/compares", Dict(:add => [other], :subtraction => other))
    _workload_call(send, "PATCH", "/api/samples/$sample/compares", Dict(:remove => [other]))
    _workload_call(send, "GET", "/api/experiments/$experiment/labels")
    _workload_call(send, "POST", "/api/samples/$sample/labels", Dict(:label => "precompile"))
    _workload_call(send, "POST", "/api/samples/$sample/labels/suggested", Dict(:label => "precompile"))
    _workload_call(send, "DELETE", "/api/samples/$sample/labels?label=precompile")
    _workload_call(send, "GET", "/api/label-colors")
    _workload_call(send, "PUT", "/api/label-colors", Dict(:label => "precompile", :color => "teal"))
    _workload_call(send, "PUT", "/api/label-colors", Dict(:label => "precompile", :color => nothing))
    _workload_call(send, "GET", "/api/exposures/$exposure/trace")
    _workload_call(send, "GET", "/api/exposures/$exposure/image?thumb=1")
    _workload_call(send, "GET", "/api/exposures/$exposure/image")
    peaks = _workload_call(send, "GET", "/api/exposures/$exposure/peaks")
    indices = _workload_call(send, "GET", "/api/exposures/$exposure/indices")
    _workload_call(send, "GET", "/api/exposures/$exposure/assignment")
    for part in ("trace", "peaks", "indices", "assignment")
        _workload_call(send, "GET", "/api/exposures/$mixed/$part")
    end

    added = _workload_call(send, "POST", "/api/exposures/$exposure/peaks", Dict(:q => 0.3))
    _workload_call(send, "DELETE", "/api/peaks/manual/$(added.id)")
    at = findfirst(p -> p.source == "auto", peaks)
    at === nothing && error("precompile workload: the synthetic trace gave no auto peaks")
    auto = peaks[at]
    _workload_call(send, "PATCH", "/api/peaks/auto/$(auto.id)", Dict(:excluded => true))
    _workload_call(send, "PATCH", "/api/peaks/auto/$(auto.id)", Dict(:excluded => false))

    _workload_call(send, "PATCH", "/api/exposures/$exposure/status", Dict(:status => "rejected"))
    _workload_call(send, "PATCH", "/api/exposures/$exposure/status", Dict(:status => nothing))
    _workload_call(send, "PUT", "/api/samples/$sample/selection", Dict(:exposure_ids => Int[]))
    _workload_call(send, "PUT", "/api/samples/$sample/selection", Dict(:exposure_ids => [exposure]))

    indices = _workload_call(send, "GET", "/api/exposures/$exposure/indices")
    if !isempty(indices)
        index_id = first(indices).id
        _workload_call(send, "POST", "/api/exposures/$exposure/assignment/members", Dict(:index_id => index_id))
        _workload_call(send, "DELETE", "/api/exposures/$exposure/assignment/members/$index_id")
    end
    _workload_call(send, "GET", "/api/samples/$sample/call")
    _workload_call(send, "POST", "/api/samples/$sample/call/state", Dict(:state => "form_factor"))
    _workload_call(send, "POST", "/api/samples/$sample/call/state", Dict(:state => "indexed"))
    _workload_call(send, "POST", "/api/exposures/$exposure/custom-index",
                   Dict(:phase => "Pn3m", :basis => 0.06, :ratios => sqrt.(Float64[2, 3, 4])))
    for phase in _WORKLOAD_PHASES
        _workload_call(send, "POST", "/api/exposures/$exposure/custom-index",
                       Dict(:phase => string(nameof(phase)), :basis => 0.06))
    end
    autos = [p for p in peaks if p.source == "auto"]
    if length(autos) >= 2
        spec = _workload_call(send, "POST", "/api/exposures/$exposure/speculative",
                              Dict(:phase => "Pn3m", :anchor_peak_id => autos[1].id, :anchor_ratio => 1,
                                   :additional => [Dict(:ratio_position => 2, :peak_id => autos[2].id)]))
        _workload_call(send, "DELETE", "/api/indices/$(spec.id)")
    end
    _workload_call(send, "PATCH", "/api/experiments/$experiment", Dict(:name => "precompile renamed", :description => "precompile"))
    _workload_call(send, "PATCH", "/api/experiments/$experiment", Dict(:name => "precompile", :description => nothing))
    for setup in get(_workload_call(send, "GET", "/api/experiments/$experiment"), :setups, [])
        _workload_call(send, "PATCH", "/api/setups/$(setup.id)", Dict(:name => "precompile setup"))
        _workload_call(send, "PATCH", "/api/setups/$(setup.id)", Dict(:name => nothing))
    end
    _workload_call(send, "GET", "/api/experiments/$experiment/export")
    _workload_call(send, "GET", "/api/experiments/$experiment/export?format=csv")
    _workload_call(send, "POST", "/api/exposures/$exposure/analyze", Dict())
end

function _workload_peak_model!(db::SQLite.DB, experiment::Int, sample::Int, image::AbstractString, analysis_dir::AbstractString)
    _workload_trace(joinpath(analysis_dir, replace(DEFAULT_INTEGRATION_PATTERN, "{name}" => "frame_model")))
    exposure = create_exposure!(db; experiment_id = experiment, sample_id = sample,
                                filename = "frame_model", image_path = image)
    analyze_exposure!(db, exposure, analysis_dir)
    recipe = HimalayaPeakModel.default_recipe(n_trees = 3, threshold = 0.3)
    rng = HimalayaPeakModel.MersenneTwister(1)
    X = [mod(0.37i + 0.11j + 0.013i * j, 1.0) for i in 1:80, j in eachindex(HimalayaPeakModel.feature_names(recipe.features))]
    fitted = HimalayaPeakModel.fit(recipe.learner, X, X[:, 1] .> 0.5; rng)
    id = insert_peak_model!(db, Dict{String, Any}("status" => "passed", "watermark" => 0,
        "counts" => Dict("curations" => 1), "recipe" => HimalayaPeakModel.to_dict(recipe),
        "recipe_id" => HimalayaPeakModel.recipe_id(recipe),
        "feature_version" => HimalayaPeakModel.feature_version(recipe.features),
        "core_version" => string(pkgversion(Himalaya)), "threshold" => recipe.threshold, "nms" => recipe.nms,
        "metrics" => Dict{String, Any}(), "fitted" => HimalayaPeakModel.to_dict(fitted)))
    promote_peak_model!(db, id) || error("precompile workload: the synthetic peak model was not promoted")
    isempty(retire_stale_peak_models!(db)) || error("precompile workload: the synthetic peak model was retired")
    stopped = peak_model_pass_stopped()
    set_peak_model_pass_stopped!(false)
    try
        run_peak_model_pass!(db)
    finally
        set_peak_model_pass_stopped!(stopped)
    end
    only(sqlrows(DBInterface.execute(db, "SELECT peak_model_version AS v FROM exposures WHERE id = ?", [exposure]))).v == id ||
        error("precompile workload: the peak model did not write the synthetic exposure's peaks")
    nothing
end

function _save_module_state()
    (db = _DB_REF[],
     pool = _READER_POOL[],
     subscribers = lock(() -> copy(SSE_SUBSCRIBERS[]), SSE_LOCK),
     op_ids = lock(() -> Set(keys(OP_LOCKS)), OP_LOCKS_MU),
     models = lock(() -> copy(PEAK_MODEL_CACHE), PEAK_MODEL_CACHE_LOCK),
     hashes = lock(() -> (copy(_PEAK_SET_HASH_CACHE), copy(_PEAK_SET_HASH_CACHE_ORDER)), _PEAK_SET_HASH_CACHE_LOCK),
     cache_bytes = lock(() -> copy(_CACHE_BYTES), _IMAGE_CACHE_LOCK),
     prp = lock(() -> copy(PRP_DETAILS), PRP_DETAILS_LOCK))
end

function _restore_module_state!(saved)
    _close_reader_pool!(_READER_POOL[])
    _DB_REF[] = saved.db
    _READER_POOL[] = saved.pool
    lock(() -> SSE_SUBSCRIBERS[] = saved.subscribers, SSE_LOCK)
    lock(() -> filter!(kv -> kv.first in saved.op_ids, OP_LOCKS), OP_LOCKS_MU)
    lock(() -> (empty!(PEAK_MODEL_CACHE); merge!(PEAK_MODEL_CACHE, saved.models)), PEAK_MODEL_CACHE_LOCK)
    lock(_PEAK_SET_HASH_CACHE_LOCK) do
        empty!(_PEAK_SET_HASH_CACHE); merge!(_PEAK_SET_HASH_CACHE, saved.hashes[1])
        empty!(_PEAK_SET_HASH_CACHE_ORDER); append!(_PEAK_SET_HASH_CACHE_ORDER, saved.hashes[2])
    end
    lock(() -> (empty!(_CACHE_BYTES); merge!(_CACHE_BYTES, saved.cache_bytes)), _IMAGE_CACHE_LOCK)
    lock(() -> (empty!(PRP_DETAILS); merge!(PRP_DETAILS, saved.prp)), PRP_DETAILS_LOCK)
    nothing
end

function run_route_workload(dir::AbstractString; record = nothing)
    analysis_dir = joinpath(dir, "analysis")
    data_dir = joinpath(dir, "data")
    mkpath(analysis_dir); mkpath(data_dir)
    _workload_trace(joinpath(analysis_dir, replace(DEFAULT_INTEGRATION_PATTERN, "{name}" => "frame")))
    _workload_trace(joinpath(analysis_dir, replace(DEFAULT_INTEGRATION_PATTERN, "{name}" => "frame_mixed")),
                    _workload_mixed_centers())
    image = joinpath(data_dir, "frame.tif")
    _workload_frame(image)
    prp = joinpath(data_dir, "frame.prp")
    _workload_prp(prp)
    saved = _save_module_state()
    db = open_db(joinpath(dir, "himalaya.db"))
    try
        experiment = init_experiment!(db; path = dir, data_dir = data_dir, analysis_dir = analysis_dir)
        sample = create_sample!(db; experiment_id = experiment, name = "precompile")
        other = create_sample!(db; experiment_id = experiment, name = "precompile_b")
        exposure = create_exposure!(db; experiment_id = experiment, sample_id = sample,
                                    filename = "frame", image_path = image, prp_path = prp,
                                    exposure_time = 1.0, transmitted_counts = 1000.0,
                                    incident_monitor = 1.5, counting_time = 1.05)
        create_exposure!(db; experiment_id = experiment, sample_id = sample,
                         filename = "frame_b", image_path = image, prp_path = prp,
                         timestamp = "2026-01-01T00:00:00", exposure_time = 1.0, horizontal_position = 1.0,
                         scan_id = 1, frame_no = 1, transmitted_counts = 1.0, incident_monitor = 1.0,
                         counting_time = 1.0)
        analyze_exposure!(db, exposure, analysis_dir)
        mixed = create_exposure!(db; experiment_id = experiment, sample_id = other,
                                 filename = "frame_mixed", image_path = image)
        analyze_exposure!(db, mixed, analysis_dir)

        Oxygen.resetstate()
        _DB_REF[] = db
        _READER_POOL[] = _new_reader_pool(db)
        register_routes!()

        inproc = _recording(_workload_inproc, record)
        _workload_call(inproc, "POST", "/api/users", Dict(:username => "precompile"))
        _workload_requests(inproc, experiment, sample, other, exposure, mixed)
        _workload_peak_model!(db, experiment, sample, image, analysis_dir)
        _gc_sweep(db, 3600)
        _boot_background_work(db, false)
        port = find_free_port()
        Oxygen.serve(; host = "127.0.0.1", port, async = true, show_banner = false,
                     docs = false, metrics = false, parallel = true,
                     middleware = [db_connection_middleware])
        try
            _wait_for_server(port)
            _workload_requests(_recording(_workload_socket(port), record), experiment, sample, other, exposure, mixed)
        finally
            Oxygen.terminate()
        end
    finally
        Oxygen.resetstate()
        _restore_module_state!(saved)
        DBInterface.close!(db)
    end
    nothing
end

@setup_workload begin
    @compile_workload begin
        Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            parse_serve_args(["--port", "8080", "--host", "127.0.0.1"])
            mktempdir(run_route_workload)
        end
    end
end
