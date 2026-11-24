using Printf, Statistics

const USAGE = """
usage: julia packages/HimalayaUI/scripts/check_transmitted_beam.jl <beamtime folder> [--outliers N]

Reads every .prp under the folder and changes nothing. For each folder of .prp files
(one per setup) it prints:
  - how many images name each I_Signal channel
  - how many have no usable I_2, so subtraction falls back to the high-q match
  - how many look beam-off (I_0 at or below zero)
  - the spread of I_2 / I_signal against the folder median. Scaling a buffer by I_2
    instead of I_signal changes its scale by the ratio of the two images' values.
  - I_2 on beam-off images, which should be near zero
  - the spread of I_2 / I_0 against the folder median, with images more than
    STEP_FRACTION away listed: a jump there is a gain or attenuator step, or a
    blocked or moved beamstop
  - I_2 per second of counting time for each counting time seen, which should agree
"""

const STEP_FRACTION = 0.2

number(line) = (m = match(r"=\s*([+-]?[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)", line); m === nothing ? nothing : parse(Float64, m[1]))

function read_monitors(path)
    channel, signal, i2, i0, counting = nothing, nothing, nothing, nothing, nothing
    for line in eachline(path)
        if startswith(line, "I_Signal (i.e. I1 in header file) is from:")
            channel = strip(split(line, ":"; limit = 2)[2])
        elseif startswith(lowercase(line), "i_signal=")
            signal = number(line)
        elseif startswith(line, "I_2(")
            i2 = number(line)
        elseif startswith(line, "I_0(")
            i0 = number(line)
        elseif startswith(line, "Counting time=")
            counting = number(line)
        end
    end
    (; channel, signal, i2, i0, counting)
end

usable(x) = x !== nothing && isfinite(x) && x > 0

function prp_folders(root)
    found = Dict{String,Vector{String}}()
    for (dir, _, files) in walkdir(root)
        prps = [joinpath(dir, f) for f in files if endswith(f, ".prp")]
        isempty(prps) || (found[dir] = sort(prps))
    end
    found
end

spread(xs) = quantile(xs, [0.0, 0.05, 0.5, 0.95, 1.0])

function report_beam_off(io, rows)
    i2 = [m.i2 for (_, m) in rows if m.i2 !== nothing]
    isempty(i2) && return
    @printf(io, "   I_2 on beam-off images: min %.1f  median %.1f  max %.1f\n", minimum(i2), median(i2), maximum(i2))
end

function report_incident_ratio(io, lit, outliers)
    ratio = [(n, m.i2 / m.i0) for (n, m) in lit if usable(m.i2) && usable(m.i0)]
    isempty(ratio) && return
    mid = median(last.(ratio))
    rel = [(n, r / mid) for (n, r) in ratio]
    @printf(io, "   I_2 / I_0 vs folder median: min %.3f  p5 %.3f  median %.3f  p95 %.3f  max %.3f\n", spread(last.(rel))...)
    steps = sort([x for x in rel if abs(x[2] - 1) > STEP_FRACTION]; by = x -> abs(log(x[2])), rev = true)
    @printf(io, "   more than %d%% off (gain, attenuator or beamstop step?): %d\n", round(Int, 100 * STEP_FRACTION), length(steps))
    for (n, r) in first(steps, outliers)
        @printf(io, "      %-48s %.3f\n", n, r)
    end
end

function report_counting_time(io, lit)
    by_time = Dict{Float64,Vector{Float64}}()
    for (_, m) in lit
        (usable(m.i2) && usable(m.counting)) || continue
        push!(get!(by_time, m.counting, Float64[]), m.i2 / m.counting)
    end
    isempty(by_time) && return
    println(io, "   I_2 per counting second, by counting time:")
    for t in sort(collect(keys(by_time)))
        rates = by_time[t]
        @printf(io, "      %8.2f s  × %-5d median %.4g\n", t, length(rates), median(rates))
    end
end

function report_folder(io, root, dir, paths, outliers)
    rows = [(basename(p), read_monitors(p)) for p in paths]
    channels = Dict{String,Int}()
    for (_, m) in rows
        key = something(m.channel, "(no I_Signal line)")
        channels[key] = get(channels, key, 0) + 1
    end
    beam_off = count(((_, m),) -> m.i0 !== nothing && m.i0 <= 0, rows)
    lit = [(n, m) for (n, m) in rows if !(m.i0 !== nothing && m.i0 <= 0)]
    fallback = [n for (n, m) in lit if !usable(m.i2)]
    paired = [(n, m.i2 / m.signal) for (n, m) in lit if usable(m.i2) && usable(m.signal)]

    println(io, "\n== ", relpath(dir, root), "  (", length(rows), " .prp)")
    println(io, "   I_Signal channel: ", join(["$k × $v" for (k, v) in sort(collect(channels))], ", "))
    println(io, "   beam-off (I_0 ≤ 0): ", beam_off)
    println(io, "   no usable I_2, falls back to high-q match: ", length(fallback), " of ", length(lit))
    for n in first(fallback, outliers)
        println(io, "      ", n)
    end
    report_beam_off(io, [(n, m) for (n, m) in rows if m.i0 !== nothing && m.i0 <= 0])
    report_incident_ratio(io, lit, outliers)
    report_counting_time(io, lit)
    isempty(paired) && return
    mid = median(last.(paired))
    rel = [(n, r / mid) for (n, r) in paired]
    q = spread(last.(rel))
    @printf(io, "   I_2 / I_signal vs folder median: min %.3f  p5 %.3f  median %.3f  p95 %.3f  max %.3f\n", q...)
    @printf(io, "   so a buffer scale would move by up to %+.1f%% / %+.1f%% (p5 / p95 pair against the median)\n",
            100 * (q[2] - 1), 100 * (q[4] - 1))
    sorted = sort(rel; by = r -> abs(log(r[2])), rev = true)
    println(io, "   largest departures:")
    for (n, r) in first(sorted, outliers)
        @printf(io, "      %-48s %.3f\n", n, r)
    end
end

function main(args)
    isempty(args) && (print(USAGE); return 1)
    root = abspath(args[1])
    isdir(root) || (println(stderr, "not a folder: ", root); return 1)
    outliers = 5
    i = findfirst(==("--outliers"), args)
    i === nothing || (outliers = parse(Int, args[i + 1]))
    folders = prp_folders(root)
    isempty(folders) && (println("no .prp files under ", root); return 1)
    for dir in sort(collect(keys(folders)))
        report_folder(stdout, root, dir, folders[dir], outliers)
    end
    0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
