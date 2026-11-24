const AUTOSAMPLER_DIR = "autosampler"
const AUTOSAMPLER_TRUTH_FILE = "TRUTH_autosampler.csv"
const AUTOSAMPLER_TRUTH_HEADER = "sample,stem,run,load,well,kind,default_background,frames_total,frames_rejected,flags,repeat_of,radius_A,amplitude,minimum_q,note"
const AS_FRAMES = 10
const AS_EXPOSURE_S = 1
const AS_HORIZ_MM = 61.2
const AS_FIRST_SCAN = 101
const AS_CALIBRATION_SCAN = 100
const AS_T0 = DateTime(2026, 4, 27, 9, 30, 0)
const AS_LOG_AGE = Hour(1)
const AS_COLLAPSED_FRACTION = 0.04
const AS_SPHERE_FIRST_MINIMUM = 4.4934

Base.@kwdef struct AsShot
    well::Char
    kind::Char = 'S'
    radius::Float64 = 0.0
    amplitude::Float64 = 0.0
    rejected::Int = 0
    collapsed::Bool = false
end

const AS_RUNS = [
    (label = "lysA", plate = "01", queue = true, shots = [
        AsShot(well = 'A', kind = 'B'),
        AsShot(well = 'B', radius = 28.0, amplitude = 6.0),
        AsShot(well = 'C', radius = 34.0, amplitude = 8.0, rejected = 1),
        AsShot(well = 'D', radius = 40.0, amplitude = 5.0, rejected = 4),
        AsShot(well = 'E', radius = 46.0, amplitude = 9.0),
        AsShot(well = 'F', radius = 30.0, amplitude = 7.0, collapsed = true),
        AsShot(well = 'H', kind = 'B'),
    ]),
    (label = "micB", plate = "02", queue = false, shots = [
        AsShot(well = 'A', kind = 'B', rejected = 1),
        AsShot(well = 'B', radius = 52.0, amplitude = 4.0),
        AsShot(well = 'C', radius = 56.0, amplitude = 6.0, rejected = 9),
        AsShot(well = 'D', radius = 60.0, amplitude = 5.0, rejected = 1),
        AsShot(well = 'E', radius = 64.0, amplitude = 8.0),
    ]),
    (label = "lysA", plate = "01", queue = true, shots = [
        AsShot(well = 'A', kind = 'B'),
        AsShot(well = 'C', radius = 34.0, amplitude = 7.6),
        AsShot(well = 'G', radius = 38.0, amplitude = 4.0, rejected = 1),
        AsShot(well = 'H', kind = 'B'),
    ]),
]

as_scan(n) = lpad(n, 4, '0')
as_stem(label, plate, shot::AsShot, scan) = "$(label)_$(plate)$(shot.well)_$(shot.kind)$(as_scan(scan))"
as_name(label, plate, shot::AsShot, scan) = "$label $plate$(shot.well) ($(shot.kind)$(as_scan(scan)))"

as_dirs(root) = (data = joinpath(root, "data"),
                 analysis = joinpath(root, "analysis"),
                 tot = joinpath(root, "analysis", "automatic_analysis", "tot_files"),
                 logs = joinpath(root, "analysis", "automatic_analysis", "sastool"))

function as_prepare(root)
    rm(root; force = true, recursive = true)
    dirs = as_dirs(root)
    foreach(mkpath, (dirs.data, dirs.tot, dirs.logs))
    dirs
end

function polydisperse_sphere_form_factor(q, radius)
    zs = (-1.5, -0.75, 0.0, 0.75, 1.5)
    ws = exp.(-0.5 .* zs .^ 2)
    sum(w * begin
        x = q * radius * (1 + 0.07 * z)
        x < 1e-6 ? 1.0 : (3 * (sin(x) - x * cos(x)) / x^3)^2
    end for (w, z) in zip(ws, zs)) / sum(ws)
end

function as_buffer_model(rise_scale)
    @. SOLVENT_FLOOR * (1 + RISE * rise_scale * Q^-RISE_EXP)
end

function as_model(buffer, shot::AsShot)
    shot.kind == 'B' && return copy(buffer)
    buffer .+ shot.amplitude * SOLVENT_FLOOR .* polydisperse_sphere_form_factor.(Q, shot.radius)
end

function as_rejected_frames(total, rejected)
    rejected == 0 && return Int[]
    2 .+ floor.(Int, (0:rejected-1) .* (total - 1) ./ rejected)
end

function as_frame_line(rng, i, rejected::Bool)
    variance = rejected ? 140 + 60 * rand(rng) : 20 + 70 * rand(rng)
    @sprintf("frame #%d: %.4f  %.2f  %.1f   -nan%s\n", i, 0.0004 + 0.0008 * rand(rng),
             -10 + 0.4 * randn(rng), variance, rejected ? "  Frame rejected!!" : "")
end

function as_frame_section(rng, total, rejected)
    bad = Set(as_rejected_frames(total, rejected))
    join(as_frame_line(rng, i, i in bad) for i in 1:total)
end

const AS_BUFFER_HEADER = "Variance of frames relative to the first frame for buffer files:\n"
const AS_SAMPLE_HEADER = "Variance of frames relative to the first frame for sample files:\n"

as_first_buffer_log(rng, total, rejected) = AS_BUFFER_HEADER * as_frame_section(rng, total, rejected)

function as_sample_log(rng, total, rejected)
    AS_BUFFER_HEADER * as_frame_section(rng, 3, 0) *
        "Average variance for buffer:\n100\n\nVariance used for rejection: 100\nFactor used for rejection: 2\n" *
        AS_SAMPLE_HEADER * as_frame_section(rng, total, rejected) *
        "Buffer frames:$total Rejected:0\nSample frames:$total Rejected:$rejected\n"
end

function as_write_log(path, text)
    write(path, text)
    run(`touch -t $(Dates.format(now() - AS_LOG_AGE, "yyyymmddHHMM.SS")) $path`)
end

function as_write_trace(path, I, sig)
    open(path, "w") do io
        for j in eachindex(Q)
            @printf(io, "%.6e %.6e %.6e\n", Q[j], I[j], sig[j])
        end
    end
end

function as_write_frames!(rng_beam, dirs, stem, t, i0, transmission, i2_transmission, emit_first, images)
    for frame in 1:AS_FRAMES
        name = "$(stem)_0_$(lpad(frame, 3, '0'))"
        write_prp(joinpath(dirs.data, name * ".prp"), Dates.format(t + Second(2 * frame), "dd u yyyy HH:MM:SS");
                  exposure_time = AS_EXPOSURE_S, horiz = AS_HORIZ_MM, i0, transmission = i2_transmission, rng = rng_beam)
        tif = joinpath(dirs.data, name * ".tif")
        if frame == 1
            emit_first(tif)
        else
            images ? hardlink(joinpath(dirs.data, "$(stem)_0_001.tif"), tif) : touch(tif)
        end
    end
end

function as_write_acquisition!(rngs, dirs, run, shot, scan, t, buffer, gen)
    stem = as_stem(run.label, run.plate, shot, scan)
    first_frame = stem * "_0_001"
    counting = counting_time(AS_EXPOSURE_S)
    f = counting / 15.8
    i0 = (1.0 - 0.18 * mod(Dates.value(t - AS_T0) / 60_000, 120) / 120) * (1 + 0.01 * randn(rngs.beam))
    transmission = 0.85 * (1 + 0.01 * randn(rngs.beam))
    i2_transmission = shot.collapsed ? AS_COLLAPSED_FRACTION * transmission : transmission
    kept = AS_FRAMES - shot.rejected
    model = as_model(buffer, shot)
    sig = @. sqrt(max(model, 1e-6) / (PX_PER_TRACE * f * kept * NPIX))
    I = @. max(model + randn(rngs.trace) * sig, 1e-6)
    flux = i0 * transmission / 0.86
    divide = f * GRAND_BASE / i2_reading(first_frame, i0, i2_transmission; exposure_time = AS_EXPOSURE_S)
    as_write_trace(joinpath(dirs.tot, stem * "_tot.dat"), I .* (flux * divide), sig .* (sqrt(flux) * divide))
    emit_first(tif) = gen.images ?
        write_tif(tif, make_frame(rngs.img, model, gen.qm, gen.dx, gen.dy, gen.gaps, f * flux)) : touch(tif)
    as_write_frames!(rngs.beam, dirs, stem, t, i0, transmission, i2_transmission, emit_first, gen.images)
    stem
end

function as_write_calibration!(rngs, dirs, buffer, gen)
    stem = "agbe_S$(as_scan(AS_CALIBRATION_SCAN))"
    write_prp(joinpath(dirs.data, stem * "_0_001.prp"), Dates.format(AS_T0 - Minute(20), "dd u yyyy HH:MM:SS");
              exposure_time = AS_EXPOSURE_S, horiz = AS_HORIZ_MM, rng = rngs.beam)
    tif = joinpath(dirs.data, stem * "_0_001.tif")
    gen.images ? write_tif(tif, make_frame(rngs.img, buffer, gen.qm, gen.dx, gen.dy, gen.gaps, 0.2)) : touch(tif)
    stem
end

function as_write_log_for(rngs, dirs, stem, first_buffer, shot)
    path = joinpath(dirs.logs, stem * "_0_001.log")
    as_write_log(path, first_buffer ? as_first_buffer_log(rngs.log, AS_FRAMES, shot.rejected) :
                                      as_sample_log(rngs.log, AS_FRAMES, shot.rejected))
end

function as_flags(run, shot, first_buffer, repeated)
    flags = String[]
    first_buffer && push!(flags, "first_buffer")
    shot.rejected >= 2 && shot.rejected < AS_FRAMES - 1 && push!(flags, "several_rejected")
    shot.rejected == AS_FRAMES - 1 && push!(flags, "one_frame_left")
    shot.collapsed && push!(flags, "i2_collapsed")
    last(run.shots).kind == 'B' && length(run.shots) > 1 || push!(flags, "no_closing_buffer")
    repeated && push!(flags, "repeated_well")
    flags
end

function as_truth_row(name, stem, run, load, shot, default_bg, flags, repeat_of)
    kind = shot.kind == 'B' ? "buffer" : "sample"
    minimum_q = shot.kind == 'B' ? "" : @sprintf("%.4f", AS_SPHERE_FIRST_MINIMUM / shot.radius)
    note = shot.collapsed ? "I_2 collapsed to a few percent of normal; the beamline's division inflates the trace" :
           shot.rejected == AS_FRAMES - 1 ? "all but one frame rejected" :
           shot.rejected >= 2 ? "several frames rejected" : ""
    join([csvq(name), stem, run.label, string(load), string(run.plate, shot.well), kind, csvq(default_bg),
          string(AS_FRAMES), string(shot.rejected), join(flags, " "), csvq(repeat_of),
          shot.kind == 'B' ? "" : string(shot.radius), shot.kind == 'B' ? "" : string(shot.amplitude),
          minimum_q, csvq(note)], ",")
end

function as_write_queue_file(dirs, run, scan)
    open(joinpath(dirs.data, "AT_SolSAXS_$(run.label)_$(as_scan(scan)).txt"), "w") do io
        println(io, "queue for series $(run.label), first scan $(as_scan(scan))")
        println(io, "load 30 uL, rinse 2 x 60 uL")
    end
end

function write_autosampler!(out, seed, images)
    root = joinpath(out, AUTOSAMPLER_DIR)
    dirs = as_prepare(root)
    write_setup_info(joinpath(dirs.analysis, "setup_info_20260301_120000.txt");
                     calibration = "../data/agbe_S$(as_scan(AS_CALIBRATION_SCAN))_0_001.tif")
    write_integ_mpp(joinpath(dirs.analysis, "integ.mpp"); normalised = true)
    rngs = (trace = Xoshiro(hash((seed, :autosampler, :trace))),
            img = Xoshiro(hash((seed, :autosampler, :image))),
            beam = Xoshiro(hash((seed, :autosampler, :beam))),
            log = Xoshiro(hash((seed, :autosampler, :log))),
            run = Xoshiro(hash((seed, :autosampler, :run))))
    qm, dx, dy = images ? q_map() : (nothing, nothing, nothing)
    gen = (images = images, qm = qm, dx = dx, dy = dy, gaps = gap_mask())
    truth = String[AUTOSAMPLER_TRUTH_HEADER]
    seen = Dict{Tuple{String,String},String}()
    t = AS_T0
    scan = AS_FIRST_SCAN
    for run in AS_RUNS
        run.queue && as_write_queue_file(dirs, run, scan)
        buffer = as_buffer_model(exp(0.2 * randn(rngs.run)))
        load = scan
        first_buffer_name = ""
        for (k, shot) in enumerate(run.shots)
            t += Second(rand(rngs.run, 75:105))
            stem = as_write_acquisition!(rngs, dirs, run, shot, scan, t, buffer, gen)
            name = as_name(run.label, run.plate, shot, scan)
            first_buffer = shot.kind == 'B' && isempty(first_buffer_name)
            first_buffer && (first_buffer_name = name)
            as_write_log_for(rngs, dirs, stem, first_buffer, shot)
            well = (run.label, string(run.plate, shot.well))
            repeat_of = shot.kind == 'S' ? get(seen, well, "") : ""
            shot.kind == 'S' && !haskey(seen, well) && (seen[well] = name)
            default_bg = shot.kind == 'B' ? "" : first_buffer_name
            push!(truth, as_truth_row(name, stem, run, load, shot, default_bg,
                                      as_flags(run, shot, first_buffer, !isempty(repeat_of)), repeat_of))
            scan += 1
        end
        t += Minute(run.queue ? 25 : 4)
    end
    stem = as_write_calibration!(rngs, dirs, as_buffer_model(1.0), gen)
    push!(truth, join([stem, stem, "", "0", "", "calibration", "", "", "", "calibration", "", "", "", "",
                       csvq("calibration shot, no averaged trace")], ","))
    write(joinpath(out, AUTOSAMPLER_TRUTH_FILE), join(truth, "\n") * "\n")
    (acquisitions = length(truth) - 2, runs = length(AS_RUNS), calibration = 1)
end
