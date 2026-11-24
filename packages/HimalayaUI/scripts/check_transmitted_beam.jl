using Printf, Statistics

const USAGE = """
usage: julia packages/HimalayaUI/scripts/check_transmitted_beam.jl <beamtime folder> [--outliers N]
       julia packages/HimalayaUI/scripts/check_transmitted_beam.jl <beamtime folder> --pair <sample stem> <buffer stem>

Reads .prp, _tot.dat and integ.mpp files under the folder and changes nothing.
It first lists every integ.mpp with its -i line: `-i yes i2 <base>` means SASTool
divided each trace by its I_2, so a buffer scale within that setup is 1. Then, for
each folder of .prp files (one per setup), it prints:
  - how many images name each I_Signal channel
  - how many have no usable I_2, so subtraction falls back to the high-q match
  - how many look beam-off (I_0 at or below zero)
  - the spread of I_2 / I_signal against the folder median. On raw-count traces,
    scaling a buffer by I_2 instead of I_signal changes its scale by this ratio.
  - I_2 on beam-off images, which should be near zero
  - the spread of I_2 / I_0 against the folder median, with images more than
    STEP_FRACTION away listed: a jump there is a gain or attenuator step, or a
    blocked or moved beamstop
  - I_2 per second of counting time for each counting time seen, which should agree
  - the log-log slope of each trace's mid-q σ²/I against its I_2: near −1 when the
    traces are divided by I_2, near 0 when they are raw counts
  - the spread of I_signal against the folder median: on I_2-divided traces the
    released app scales a buffer by the I_signal ratio, so this is its error

--pair prints, for one sample and its buffer, the I_2 ratio and the ratio of their
high-q levels at scale 1. On I_2-divided traces the second should be near 1 even
when the first is not.
"""

const HIGH_Q_FRACTION = 0.15

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

function integ_lines(root)
    found = Tuple{String,String}[]
    for (dir, _, files) in walkdir(root)
        "integ.mpp" in files || continue
        path = joinpath(dir, "integ.mpp")
        line = something(findfirst(l -> startswith(strip(l), "-i"), readlines(path)), 0)
        push!(found, (relpath(path, root), line == 0 ? "(no -i line: raw counts)" : strip(readlines(path)[line])))
    end
    found
end

function tot_files(root)
    found = Dict{String,String}()
    for (dir, _, files) in walkdir(root), f in files
        endswith(f, "_tot.dat") && (found[f[1:end-length("_tot.dat")]] = joinpath(dir, f))
    end
    found
end

tot_key(prp_name) = replace(prp_name[1:end-length(".prp")], r"_\d+_\d+$" => "")

function high_q_level(path)
    rows = [parse.(Float64, split(l)[1:2]) for l in eachline(path) if length(split(l)) >= 2 && tryparse(Float64, split(l)[1]) !== nothing]
    isempty(rows) && return nothing
    qs = first.(rows)
    from = maximum(qs) - (maximum(qs) - minimum(qs)) * HIGH_Q_FRACTION
    levels = [r[2] for r in rows if r[1] >= from && r[2] > 0]
    isempty(levels) ? nothing : median(levels)
end

slope(xs, ys) = sum(xs .* ys) / sum(abs2, xs)

# physics: at mid q the error column is the Poisson σ of the ring average, so σ²/I is
# physics: one over the counts behind it. Raw counts keep σ²/I the same whatever I_2
# physics: reads; dividing a trace by I_2 also divides σ²/I by I_2, a log-log slope of −1.
function mid_q_variance_ratio(path)
    rows = [parse.(Float64, split(l)[1:3]) for l in eachline(path)
            if length(split(l)) >= 3 && all(t -> tryparse(Float64, t) !== nothing, split(l)[1:3])]
    length(rows) < 10 && return nothing
    mid = rows[div(length(rows), 3):div(2 * length(rows), 3)]
    ratios = [r[3]^2 / r[2] for r in mid if r[2] > 0 && r[3] > 0]
    isempty(ratios) ? nothing : median(ratios)
end

function report_variance_slope(io, lit, tots)
    pts = Tuple{Float64,Float64}[]
    for (n, m) in lit
        usable(m.i2) || continue
        path = get(tots, tot_key(n), nothing)
        path === nothing && continue
        r = mid_q_variance_ratio(path)
        r === nothing || push!(pts, (log(m.i2), log(r)))
    end
    xs = first.(pts) .- mean(first.(pts)); ys = last.(pts) .- mean(last.(pts))
    (length(pts) < 3 || sum(abs2, xs) == 0) && return println(io, "   σ²/I vs I_2: too few traces found (", length(pts), ")")
    @printf(io, "   mid-q σ²/I vs I_2, log-log slope over %d traces: %.2f (≈ −1 divided by I_2, ≈ 0 raw counts)\n",
            length(pts), slope(xs, ys))
end

function report_signal_spread(io, lit)
    signal = [m.signal for (_, m) in lit if usable(m.signal)]
    isempty(signal) && return
    q = spread(signal ./ median(signal))
    @printf(io, "   I_signal vs folder median: p5 %.3f  p95 %.3f, so on I_2-divided traces today's buffer scale is off by %+.1f%% / %+.1f%%\n",
            q[2], q[4], 100 * (q[2] - 1), 100 * (q[4] - 1))
end

function report_pair(io, root, a, b)
    prps = Dict(basename(p)[1:end-4] => p for ps in values(prp_folders(root)) for p in ps)
    tots = tot_files(root)
    for stem in (a, b)
        haskey(prps, stem) || (println(io, "no .prp named ", stem); return 1)
        haskey(tots, tot_key(stem * ".prp")) || (println(io, "no _tot.dat for ", stem); return 1)
    end
    ma, mb = read_monitors(prps[a]), read_monitors(prps[b])
    la, lb = high_q_level(tots[tot_key(a * ".prp")]), high_q_level(tots[tot_key(b * ".prp")])
    @printf(io, "I_2 ratio sample / buffer:            %.4f\n", ma.i2 / mb.i2)
    @printf(io, "I_signal ratio sample / buffer:       %.4f\n", ma.signal / mb.signal)
    @printf(io, "high-q level ratio at scale 1:        %.4f\n", la / lb)
    0
end

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

function report_folder(io, root, dir, paths, outliers, tots)
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
    report_variance_slope(io, lit, tots)
    report_signal_spread(io, lit)
    isempty(paired) && return
    mid = median(last.(paired))
    rel = [(n, r / mid) for (n, r) in paired]
    q = spread(last.(rel))
    @printf(io, "   I_2 / I_signal vs folder median: min %.3f  p5 %.3f  median %.3f  p95 %.3f  max %.3f\n", q...)
    @printf(io, "   on raw-count traces, scaling by I_2 instead of I_signal moves a buffer scale by up to %+.1f%% / %+.1f%%\n",
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
    p = findfirst(==("--pair"), args)
    if p !== nothing
        p + 2 <= length(args) || (print(USAGE); return 1)
        return report_pair(stdout, root, args[p + 1], args[p + 2])
    end
    outliers = 5
    i = findfirst(==("--outliers"), args)
    if i !== nothing
        n = i < length(args) ? tryparse(Int, args[i + 1]) : nothing
        n === nothing && (print(USAGE); return 1)
        outliers = n
    end
    folders = prp_folders(root)
    isempty(folders) && (println("no .prp files under ", root); return 1)
    integs = integ_lines(root)
    println(isempty(integs) ? "no integ.mpp found: the app will match at high q" : "integ.mpp files:")
    for (path, line) in integs
        println("   ", path, ": ", line)
    end
    tots = tot_files(root)
    for dir in sort(collect(keys(folders)))
        report_folder(stdout, root, dir, folders[dir], outliers, tots)
    end
    0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
