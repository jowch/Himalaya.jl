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
"""

number(line) = (m = match(r"=\s*([+-]?[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)", line); m === nothing ? nothing : parse(Float64, m[1]))

function read_monitors(path)
    channel, signal, i2, i0 = nothing, nothing, nothing, nothing
    for line in eachline(path)
        if startswith(line, "I_Signal (i.e. I1 in header file) is from:")
            channel = strip(split(line, ":"; limit = 2)[2])
        elseif startswith(lowercase(line), "i_signal=")
            signal = number(line)
        elseif startswith(line, "I_2(")
            i2 = number(line)
        elseif startswith(line, "I_0(")
            i0 = number(line)
        end
    end
    (; channel, signal, i2, i0)
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
    isempty(paired) && return
    mid = median(last.(paired))
    rel = [(n, r / mid) for (n, r) in paired]
    q = quantile(last.(rel), [0.0, 0.05, 0.5, 0.95, 1.0])
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
