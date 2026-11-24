const _WORKLOAD_AUTOSAMPLER_SHOTS = (
    (stem = "wat_01A_B0101", i2 = 3.0e6, log = :first_buffer, rejected = 1),
    (stem = "wat_01B_S0102", i2 = 3.0e6, log = :sample, rejected = 2),
    (stem = "wat_01C_S0103", i2 = 1.2e5, log = :sample, rejected = 0),
    (stem = "wat_01D_B0104", i2 = 3.0e6, log = :sample, rejected = 0),
)

const _WORKLOAD_FRAMES_PER_SHOT = 10

function _workload_shot_prp(path::AbstractString, minute::Integer, i2::Real)
    write(path, """
        Image file name: frame.tif
        Time this file was written: 27 Apr 2026 09:$(lpad(minute, 2, '0')):00
        I_2(daqCounter1.1)=$i2 (corrected for dark counts)
        I_0(daqIon3.0)=1.5
        Exposure time=1
        Counting time=1.8
        Beam energy=9000.0 eV
        Pipe length=1700 mm
        Horizontal position=61.2 mm
        Detector=Pilatus 1M
        """)
end

function _workload_frame_log(kind::Symbol, rejected::Integer)
    lines(n) = join(("frame #$i: 0.0010  -10.00  100.0   -nan" * (i <= rejected ? "  Frame rejected!!" : "") * "\n"
                     for i in 1:n))
    buffer_header = "Variance of frames relative to the first frame for buffer files:\n"
    kind === :first_buffer && return buffer_header * lines(_WORKLOAD_FRAMES_PER_SHOT)
    buffer_header * lines(3) * "Average variance for buffer:\n100\n\n" *
        "Variance of frames relative to the first frame for sample files:\n" * lines(_WORKLOAD_FRAMES_PER_SHOT) *
        "Buffer frames:$_WORKLOAD_FRAMES_PER_SHOT Rejected:0\n" *
        "Sample frames:$_WORKLOAD_FRAMES_PER_SHOT Rejected:$rejected\n"
end

function _workload_capillary_setup(root::AbstractString)
    data = mkpath(joinpath(root, "data"))
    tot = mkpath(joinpath(root, "analysis", "automatic_analysis", "tot_files"))
    _workload_frame(joinpath(data, "cap_S0001_0_001.tif"))
    _workload_shot_prp(joinpath(data, "cap_S0001_0_001.prp"), 0, 3.0e6)
    _workload_trace(joinpath(tot, "cap_S0001_tot.dat"))
end

function _workload_autosampler_setup(root::AbstractString)
    data = mkpath(joinpath(root, "data"))
    analysis = joinpath(root, "analysis")
    tot = mkpath(joinpath(analysis, "automatic_analysis", "tot_files"))
    logs = mkpath(joinpath(analysis, "automatic_analysis", "sastool"))
    write(joinpath(analysis, "integ.mpp"), "-o automatic_analysis\n-i yes i2 2.5e7\n-s 1\n")
    write(joinpath(data, "AT_SolSAXS_wat_0101.txt"), "queue\n")
    _workload_frame(joinpath(data, "agbe_S0100_0_001.tif"))
    _workload_shot_prp(joinpath(data, "agbe_S0100_0_001.prp"), 1, 3.0e6)
    for (k, shot) in enumerate(_WORKLOAD_AUTOSAMPLER_SHOTS)
        _workload_frame(joinpath(data, "$(shot.stem)_0_001.tif"))
        _workload_shot_prp(joinpath(data, "$(shot.stem)_0_001.prp"), 1 + k, shot.i2)
        _workload_trace(joinpath(tot, "$(shot.stem)_tot.dat"), Float64[])
        log = joinpath(logs, "$(shot.stem)_0_001.log")
        write(log, _workload_frame_log(shot.log, shot.rejected))
        run(`touch -t 202601010000 $log`)
    end
end

function _workload_beamtime_tree(root::AbstractString)
    _workload_capillary_setup(joinpath(root, "capillary"))
    _workload_autosampler_setup(joinpath(root, "autosampler"))
    root
end

function _workload_refused(send, method::AbstractString, target::AbstractString, body = nothing)
    payload = body === nothing ? UInt8[] : Vector{UInt8}(JSON3.write(body))
    res = send(method, target, payload)
    res.status == 422 || error("precompile workload: $method $target returned $(res.status), expected 422")
    nothing
end

function _workload_wait_ingested(send, experiment::Integer)
    deadline = time() + 120
    while time() < deadline
        status = _workload_call(send, "GET", "/api/experiments/$experiment").ingest_status
        status == "complete" && !ingest_running(Int(experiment)) && return nothing
        status == "failed" && error("precompile workload: the beamtime ingest failed")
        sleep(0.05)
    end
    error("precompile workload: the beamtime ingest did not finish")
end

function _workload_create_beamtime(send, root::AbstractString)
    preview = _workload_call(send, "GET", "/api/fs/preview?path=$(HTTP.escapeuri(root))")
    ordered = sort(collect(preview.setups); by = s -> s.mode == "autosampler")
    setups = [Dict(:root => s.root, :data_dir => s.data_dir, :analysis_dir => s.analysis_dir,
                   :included => true, :patterns => s.patterns, :geometry => s.geometry) for s in ordered]
    created = _workload_call(send, "POST", "/api/experiments",
                             Dict(:name => "precompile beamtime", :path => root, :setups => setups))
    _workload_wait_ingested(send, created.id)
    wait_prewarm(created.id)
    Int(created.id)
end

function _workload_autosampler_requests(send, root::AbstractString)
    experiment = _workload_create_beamtime(send, root)
    _workload_call(send, "GET", "/api/experiments")
    detail = _workload_call(send, "GET", "/api/experiments/$experiment")
    mode_of = Dict(Int(s.id) => String(s.mode) for s in detail.setups)
    samples = _workload_call(send, "GET", "/api/samples?experiment_id=$experiment")
    auto = [s for s in samples if mode_of[Int(s.setup_id)] == "autosampler"]
    capillary = first(s for s in samples if mode_of[Int(s.setup_id)] == "capillary")
    for key in ("bridges", "problems", "labels", "compares")
        _workload_call(send, "GET", "/api/experiments/$experiment/$key")
    end
    first_exposure = 0
    for s in auto
        exposures = _workload_call(send, "GET", "/api/samples/$(s.id)/exposures")
        _workload_call(send, "GET", "/api/samples/$(s.id)/call")
        for e in exposures
            first_exposure == 0 && (first_exposure = Int(e.id))
            _workload_call(send, "GET", "/api/exposures/$(e.id)/image?thumb=1")
            occursin("agbe", String(e.filename)) || _workload_call(send, "GET", "/api/exposures/$(e.id)/trace")
        end
    end
    _workload_call(send, "GET", "/api/exposures/$first_exposure/image")
    background = first(s for s in auto if occursin("(S", String(s.name)))
    buffer = first(s for s in auto if occursin("(B", String(s.name)))
    _workload_call(send, "PATCH", "/api/samples/$(background.id)/compares", Dict(:subtraction => buffer.id))
    _workload_call(send, "PATCH", "/api/samples/$(background.id)/compares", Dict(:subtraction => nothing))
    _workload_call(send, "PATCH", "/api/samples/$(background.id)/compares", Dict(:subtraction => "default"))
    _workload_call(send, "GET", "/api/experiments/$experiment/compares")
    _workload_refused(send, "POST", "/api/exposures/$first_exposure/peaks", Dict(:q => 0.3))
    _workload_refused(send, "POST", "/api/samples/$(background.id)/bridge", Dict(:to_sample_id => capillary.id))
    nothing
end
