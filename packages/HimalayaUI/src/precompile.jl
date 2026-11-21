using PrecompileTools: @setup_workload, @compile_workload
using ImageCore: Gray
using ImageCore.FixedPointNumbers: Q0f31

function _workload_trace(path::AbstractString)
    q = collect(range(0.01, 0.5; length = 920))
    centers = vcat(0.05 .* Himalaya.phaseratios(Himalaya.Pn3m)[1:6],
                   0.09 .* Himalaya.phaseratios(Himalaya.Hexagonal)[1:4])
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
    n = 96
    counts = Int32[round(Int32, 2.0e3 * exp(-hypot(i - n / 2, j - n / 2) / 8)) for i in 1:n, j in 1:n]
    counts[1:4, :] .= -1
    FileIO.save(path, Gray.(reinterpret.(Q0f31, counts)))
end

function _workload_headers()
    ["Content-Type" => "application/json",
     "X-Username" => "precompile",
     "X-Client-Id" => "precompile",
     "X-Client-Op-Id" => string(rand(UInt64))]
end

function _workload_inproc(method::AbstractString, target::AbstractString, payload::Vector{UInt8})
    Oxygen.internalrequest(HTTP.Request(method, target, _workload_headers(), payload);
                           serialize = true, catch_errors = true)
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

function _workload_requests(send, experiment::Int, sample::Int, other::Int, exposure::Int)
    _workload_call(send, "GET", "/api/health")
    _workload_call(send, "GET", "/api/users")
    _workload_call(send, "GET", "/api/experiments")
    _workload_call(send, "GET", "/api/experiments/$experiment")
    _workload_call(send, "GET", "/api/experiments/$experiment/bridges")
    _workload_call(send, "GET", "/api/samples")
    _workload_call(send, "GET", "/api/samples?experiment_id=$experiment")
    _workload_call(send, "GET", "/api/samples/$sample")
    _workload_call(send, "GET", "/api/samples/$sample/exposures")
    _workload_call(send, "POST", "/api/samples/$sample/bridge", Dict(:to_sample_id => other))
    _workload_call(send, "DELETE", "/api/samples/$sample/bridge")
    _workload_call(send, "GET", "/api/exposures/$exposure")
    _workload_call(send, "GET", "/api/exposures/$exposure/trace")
    _workload_call(send, "GET", "/api/exposures/$exposure/image?thumb=1")
    _workload_call(send, "GET", "/api/exposures/$exposure/image")
    peaks = _workload_call(send, "GET", "/api/exposures/$exposure/peaks")
    indices = _workload_call(send, "GET", "/api/exposures/$exposure/indices")
    _workload_call(send, "GET", "/api/exposures/$exposure/assignment")

    added = _workload_call(send, "POST", "/api/exposures/$exposure/peaks", Dict(:q => 0.3))
    _workload_call(send, "DELETE", "/api/peaks/manual/$(added.id)")
    at = findfirst(p -> p.source == "auto", peaks)
    at === nothing && error("precompile workload: the synthetic trace gave no auto peaks")
    auto = peaks[at]
    _workload_call(send, "PATCH", "/api/peaks/auto/$(auto.id)", Dict(:excluded => true))
    _workload_call(send, "PATCH", "/api/peaks/auto/$(auto.id)", Dict(:excluded => false))

    _workload_call(send, "PATCH", "/api/exposures/$exposure/status", Dict(:status => "rejected"))
    _workload_call(send, "PATCH", "/api/exposures/$exposure/status", Dict(:status => nothing))
    _workload_call(send, "PATCH", "/api/exposures/$exposure/select", Dict())

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
    _workload_call(send, "POST", "/api/exposures/$exposure/analyze", Dict())
end

function _save_module_state()
    (db = _DB_REF[],
     subscribers = lock(() -> copy(SSE_SUBSCRIBERS[]), SSE_LOCK),
     op_ids = lock(() -> Set(keys(OP_LOCKS)), OP_LOCKS_MU),
     hashes = lock(() -> (copy(_PEAK_SET_HASH_CACHE), copy(_PEAK_SET_HASH_CACHE_ORDER)), _PEAK_SET_HASH_CACHE_LOCK))
end

function _restore_module_state!(saved)
    _DB_REF[] = saved.db
    lock(() -> SSE_SUBSCRIBERS[] = saved.subscribers, SSE_LOCK)
    lock(() -> filter!(kv -> kv.first in saved.op_ids, OP_LOCKS), OP_LOCKS_MU)
    lock(_PEAK_SET_HASH_CACHE_LOCK) do
        empty!(_PEAK_SET_HASH_CACHE); merge!(_PEAK_SET_HASH_CACHE, saved.hashes[1])
        empty!(_PEAK_SET_HASH_CACHE_ORDER); append!(_PEAK_SET_HASH_CACHE_ORDER, saved.hashes[2])
    end
    nothing
end

function run_route_workload(dir::AbstractString; record = nothing)
    analysis_dir = joinpath(dir, "analysis")
    data_dir = joinpath(dir, "data")
    mkpath(analysis_dir); mkpath(data_dir)
    _workload_trace(joinpath(analysis_dir, replace(DEFAULT_INTEGRATION_PATTERN, "{name}" => "frame")))
    image = joinpath(data_dir, "frame.tif")
    _workload_frame(image)
    saved = _save_module_state()
    db = open_db(joinpath(dir, "himalaya.db"))
    try
        experiment = init_experiment!(db; path = dir, data_dir = data_dir, analysis_dir = analysis_dir)
        sample = create_sample!(db; experiment_id = experiment, name = "precompile")
        other = create_sample!(db; experiment_id = experiment, name = "precompile_b")
        exposure = create_exposure!(db; experiment_id = experiment, sample_id = sample,
                                    filename = "frame", image_path = image)
        create_exposure!(db; experiment_id = experiment, sample_id = sample,
                         filename = "frame_b", image_path = image, prp_path = joinpath(data_dir, "frame.prp"),
                         timestamp = "2026-01-01T00:00:00", exposure_time = 1.0, horizontal_position = 1.0,
                         scan_id = 1, frame_no = 1, transmitted_counts = 1.0, incident_monitor = 1.0,
                         counting_time = 1.0)
        analyze_exposure!(db, exposure, analysis_dir)

        Oxygen.resetstate()
        bind_db!(db)
        register_routes!()

        inproc = _recording(_workload_inproc, record)
        _workload_call(inproc, "POST", "/api/users", Dict(:username => "precompile"))
        _workload_requests(inproc, experiment, sample, other, exposure)
        port = find_free_port()
        Oxygen.serve(; host = "127.0.0.1", port, async = true, show_banner = false,
                     docs = false, metrics = false, parallel = true)
        try
            _wait_for_server(port)
            _workload_requests(_recording(_workload_socket(port), record), experiment, sample, other, exposure)
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
