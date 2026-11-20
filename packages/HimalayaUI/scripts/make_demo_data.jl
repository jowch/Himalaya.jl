# Synthetic SSRL-style SAXS experiment for exercising HimalayaUI end to end.
#
# NOT REAL DATA. Every trace and detector frame is generated from the phase-ratio
# series in `Himalaya.phaseratios` plus a power-law background and counting
# noise. It exists so the ingest funnel, corpus table, loupe and focus views have
# something plausible to work on without committing lab data to the repo.
# The output is generated on demand and must never be committed (default target
# `scratch/demo-experiment/` is gitignored).
#
# Usage (from the repo root; `make demo-data` wraps this):
#   julia --project=packages/HimalayaUI packages/HimalayaUI/scripts/make_demo_data.jl \
#       [out_dir] [--seed 7] [--no-images]
#
# Layout written (what routes_fs.jl `resolve_experiment_layout` auto-detects);
# <base> = <label>_S<scan>, one scan per exposure, 1-3 scans per sample:
#   <out>/data/<base>_0_001.tif        Pilatus 1M frame, Int32 counts (-1 in module gaps)
#   <out>/data/<base>_0_001.prp        SSRL PRP sidecar (energy, pipe length, position, time)
#   <out>/analysis/setup_info_<ts>.txt beam center + AgBe-calibrated distance
#   <out>/analysis/automatic_analysis/tot_files/<base>_tot.dat   q  I  sigma
#   <out>/TRUTH.csv                    per sample: generated phases, expected auto call,
#                                      phase-tagged reflection q values
# test/test_demo_data.jl holds the regression floor on the expected auto calls.
#
# Output is deterministic for a given seed and Julia version (Xoshiro streams are
# stable within a Julia release, not guaranteed across releases).

using Himalaya: Phase, Lamellar, Hexagonal, Pn3m, Im3m, Ia3d, phaseratios
using FileIO, TiffImages, ColorTypes
using ImageCore.FixedPointNumbers: Q0f31
using Dates, Printf, Random

# ---- instrument (matches the SSRL 2026-04 defaults used in the test fixtures) ----
const ENERGY_EV        = 9000.027604502573
const WAVELENGTH_A     = 12398.42 / ENERGY_EV   # ~1.3776 Å
const PIPE_LENGTH_MM   = 1700
const MEAN_DISTANCE_MM = 1809.5                 # AgBe-calibrated sample-detector distance
const BEAM_X, BEAM_Y   = 421.409, 836.946       # detector px, origin bottom-left, y-up
const PIXEL_UM         = 172.0
const DET_W, DET_H     = 981, 1043              # Pilatus 1M: 2 x 5 modules of 487 x 195
const MOD_W, MOD_H, GAP_X, GAP_Y = 487, 195, 7, 17
const SLOT0_MM, SLOT_PITCH_MM = 58.9, 3.95      # sample-holder slot positions (real SSRL pitch)

# q grid shaped like test/data/example_tot.dat (919 points)
const Q = collect(range(5.824246e-3, 4.044763e-1; length = 919))

# Relative intensity by reflection order (lipid mesophases: the first one or two
# reflections dominate, higher orders fall off fast). The q positions themselves
# come from `phaseratios`; these only shape the heights, and orders past the end
# of a list are not drawn.
const REL = Dict{Type{<:Phase},Vector{Float64}}(
    Lamellar  => [1.0, 0.18, 0.05, 0.02, 0.008, 0.004],
    Hexagonal => [1.0, 0.22, 0.14, 0.05, 0.03, 0.012, 0.01],
    Pn3m      => [1.0, 0.75, 0.12, 0.10, 0.06, 0.05, 0.03, 0.02, 0.015],
    Im3m      => [1.0, 0.55, 0.30, 0.10, 0.06, 0.035, 0.02],
    Ia3d      => [1.0, 0.65, 0.06, 0.08, 0.04, 0.03, 0.02],
)

"q of the first allowed reflection for lattice parameter `a` (Å)."
first_peak_q(::Type{Lamellar}, a) = 2π / a                    # a = d-spacing
first_peak_q(::Type{Hexagonal}, a) = 4π / (sqrt(3) * a)       # a = rod spacing
first_peak_q(::Type{P}, a) where {P<:Phase} = 2π * first(phaseratios(P)) / a  # cubic: 2π√N/a

# (phase, lattice parameter Å, peak height relative to background at q1, peak sigma)
const Component = Tuple{DataType,Float64,Float64,Float64}

# A plausible lipid screen with known answers, mounted in two capillary holders
# (two loads, >1 h apart). Each sample is measured `nscans` times at its slot;
# every scan is one exposure (frame `_0_001`, which is the frame the ingest
# funnel's layout detection keys on).
#
# `auto` is what HimalayaUI's own pipeline (findpeaks -> indexpeaks -> auto_group)
# is expected to call, which is NOT always the generated truth: peak heights and
# lattices were tuned so each call holds across seeds (see
# test/test_demo_data.jl). Where it differs from `phases`, the note says why;
# those samples are the "promote the second hypothesis by hand" cases.
const SAMPLES = [
    # label            nscans components                                                  auto            note
    ("MO_H2O",         3, Component[(Pn3m, 102.0, 40.0, 7e-4)],                          "Pn3m",           "monoolein in excess water, clean Pn3m"),
    ("MO_H2O_37C",     2, Component[(Pn3m, 95.0, 35.0, 7e-4)],                           "Pn3m",           "Pn3m, lattice shrinks on heating"),
    ("MO_PEG5",        2, Component[(Pn3m, 108.0, 30.0, 8e-4), (Im3m, 138.0, 14.0, 8e-4)], "Im3m",          "Pn3m + Im3m coexistence (Bonnet ratio ~1.28). Pn3m √6 overlaps Im3m √10, so auto_group keeps Im3m; Pn3m is the second hypothesis to promote by hand"),
    ("MO_F127",        2, Component[(Im3m, 145.0, 30.0, 8e-4)],                          "Im3m",           "polymer-stabilised Im3m"),
    ("MO_lowH2O",      2, Component[(Ia3d, 150.0, 35.0, 7e-4)],                          "Ia3d",           "low hydration Ia3d"),
    ("DOPE_H2O",       3, Component[(Hexagonal, 74.0, 45.0, 9e-4)],                      "Hexagonal",      "DOPE inverse hexagonal"),
    ("DOPC_H2O",       2, Component[(Lamellar, 63.0, 60.0, 1.0e-3)],                     "Lamellar",       "DOPC multilamellar"),
    ("DOPC_DOPE_1to1", 2, Component[(Lamellar, 60.0, 20.0, 1.1e-3), (Hexagonal, 72.0, 22.0, 1.0e-3)], "Hexagonal+Lamellar", "La + HII coexistence"),
    ("MO_AMP_weak",    2, Component[(Pn3m, 99.0, 4.0, 9e-4)],                            "none",           "weak Pn3m near the noise floor; too weak to index automatically"),
    ("MO_AMP_10",      3, Component[(Pn3m, 96.0, 25.0, 9e-4), (Hexagonal, 53.0, 35.0, 1.0e-3)], "Hexagonal+Pn3m", "Pn3m drifting to HII under peptide"),
    ("buffer_blank",   1, Component[],                                                   "none",           "buffer only, no Bragg peaks"),
    ("vesicles_ULV",   2, Component[],                                                   "none",           "unilamellar vesicles, form factor only"),
]
const LOAD_SPLIT = 6        # samples 1..6 in holder 1, the rest in holder 2
const FIRST_SCAN = 2401

phase_name(::Type{P}) where {P<:Phase} = string(nameof(P))

"Linear interpolation on the global Q grid, clamped at the ends."
function interpq(x, ys)
    x <= Q[1] && return ys[1]
    x >= Q[end] && return ys[end]
    k = searchsortedlast(Q, x)
    t = (x - Q[k]) / (Q[k+1] - Q[k])
    ys[k] + t * (ys[k+1] - ys[k])
end

"Power-law low-q upturn + flat incoherent floor, shaped like example_tot.dat."
function background(rng, formfactor::Bool)
    amp = 0.7 + 0.6 * rand(rng)
    if formfactor
        # broad bilayer form-factor bump around 0.1-0.2 Å⁻¹ on a q^-2 upturn
        return @. amp * (3.0 * Q^-2.0 + 120 * exp(-((Q - 0.14) / 0.045)^2) + 30)
    end
    @. amp * (0.06 * Q^-3.2 + 1.2 * Q^-1.2 + 30)
end

"Noise-free model trace for one sample, plus the generated (phase, q) reflections."
function model_trace(rng, label, comps)
    bg = background(rng, startswith(label, "vesicles"))
    I = copy(bg)
    peaks = Tuple{String,Float64}[]
    for (P, a, height, sigma) in comps
        q1 = first_peak_q(P, a)
        ratios = phaseratios(P; normalize = true)
        bg1 = interpq(q1, bg)
        for (ratio, rel) in zip(ratios, REL[P])
            qk = q1 * ratio
            qk > Q[end] - 0.005 && break
            # height scales with local background so peaks stay visible at high q
            h = height * rel * bg1 * (interpq(qk, bg) / bg1)^0.35
            s = sigma * (1 + 0.4 * (qk - q1) / 0.3)      # mild broadening with q
            @. I += h * exp(-0.5 * ((Q - qk) / s)^2)
            push!(peaks, (phase_name(P), round(qk; digits = 5)))
        end
    end
    I, peaks
end

"One noisy measurement of a model trace; σ is the noise actually applied."
function measure(rng, I)
    # counting noise: relative error grows where intensity is low
    sig = @. 0.35 * (sqrt(I * 8.0) + 0.02 * I)
    @. max(I + randn(rng) * sig, 1e-3), sig
end

# ---- detector frame ----------------------------------------------------------
# Arrays are (row, col) = (DET_H, DET_W); row 1 is the TOP of the image, while
# the beam center is given y-up from the bottom-left, as in setup_info.

function gap_mask()
    m = falses(DET_H, DET_W)
    x0 = MOD_W                                     # one vertical gap (2 module columns)
    m[:, x0+1:x0+GAP_X] .= true
    for r in 1:4                                   # four horizontal gaps (5 module rows)
        y0 = r * MOD_H + (r - 1) * GAP_Y
        m[y0+1:y0+GAP_Y, :] .= true
    end
    m
end

"Per-pixel q (Å⁻¹) and beam-relative (dx, dy) in y-up pixel coordinates."
function q_map()
    qm = Matrix{Float64}(undef, DET_H, DET_W)
    dx = similar(qm); dy = similar(qm)
    for c in 1:DET_W, r in 1:DET_H
        x = c - 1; yup = DET_H - r
        dx[r, c] = x - BEAM_X; dy[r, c] = yup - BEAM_Y
        rm = hypot(dx[r, c], dy[r, c]) * PIXEL_UM * 1e-6
        tth = atan(rm / (MEAN_DISTANCE_MM * 1e-3))
        qm[r, c] = 4π * sin(tth / 2) / WAVELENGTH_A
    end
    qm, dx, dy
end

"Poisson sample; Knuth for small λ, rounded normal approximation above 30."
function poisson(rng, λ)
    λ <= 0 && return 0
    λ > 30 && return max(0, round(Int, λ + sqrt(λ) * randn(rng)))
    L = exp(-λ); k = 0; p = 1.0
    while true
        p *= rand(rng)
        p <= L && return k
        k += 1
    end
end

function make_frame(rng, Itrace, qm, dx, dy, gaps, frame_scale)
    counts = Matrix{Int32}(undef, DET_H, DET_W)
    for i in eachindex(qm)
        # radial profile -> 2D, scaled to plausible per-pixel counts for a 15 s exposure
        counts[i] = poisson(rng, interpq(qm[i], Itrace) * 0.05 * frame_scale)
        # beamstop + arm shadow (arm runs down-left from the beam centre)
        r = hypot(dx[i], dy[i])
        if r < 18 || (abs(dy[i] - 0.55 * dx[i]) < 6 && dx[i] < 0)
            counts[i] = 0
        end
    end
    for k in rand(rng, 1:length(counts), 6)        # a few hot pixels
        counts[k] = 1_000_000
    end
    counts[gaps] .= -1
    counts
end

# Pilatus frames are signed 32-bit integer TIFFs. TiffImages maps that pixel type
# to Gray{Q0f31}, and image.jl reinterprets back to Int32 counts, so writing the
# raw Int32 bits as Q0f31 produces the same on-disk format as the beamline.
write_tif(path, counts::Matrix{Int32}) = save(path, Gray.(reinterpret.(Q0f31, counts)))

function write_prp(path, ts; exposure_time = 15, horiz = 58.9)
    stem = splitext(basename(path))[1]
    open(path, "w") do io
        println(io, "Image file name: $stem.tif")
        println(io, "Time this file was written: $ts")
        println(io, "Exposure time=$exposure_time")
        println(io, "Counting time=$(exposure_time + 0.8)")
        println(io, "Beam energy=$ENERGY_EV eV")
        println(io, "Pipe length=$PIPE_LENGTH_MM mm")
        println(io, "Scan motor=motor\nScan range=delta mm")
        println(io, "Horizontal position=$horiz mm")
        println(io, "Vertical position=26.55 mm\ndispx position=9.8 mm\ndispy position=35.9 mm")
        println(io, "Detector=Pilatus 1M")
        println(io, "Detector mode=Pilatus; 0 for normal and 1 for dezingered")
        println(io, "Phi position=157\nAttenuator=0 um Al")
    end
end

function write_setup_info(path)
    open(path, "w") do io
        print(io, "File used for calibration:\n../data/agbe_1p7m_S1958_0_001.tif\n\n")
        print(io, "Beam center is at ($BEAM_X, $BEAM_Y)\n\n")
        println(io, "Calculated detector distances are:")
        for (i, d) in enumerate([1806.0, 1811.1, 1811.5])
            @printf(io, "Ring %d:\t pixel %.2f:\t %s mm\n", i, 247.82 * i, d)
        end
        println(io, "-"^58)
        println(io, "Mean distance:\t\t $MEAN_DISTANCE_MM mm")
    end
end

csvq(s) = occursin(r"[,\"]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s

function parse_args(args)
    out = joinpath(dirname(dirname(dirname(@__DIR__))), "scratch", "demo-experiment")
    seed = 7; images = true
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--seed"
            i < length(args) || error("--seed needs a value")
            seed = parse(Int, args[i+1]); i += 1
        elseif a == "--no-images"
            images = false
        elseif startswith(a, "-")
            error("unknown option $a (usage: make_demo_data.jl [out_dir] [--seed N] [--no-images])")
        else
            out = a
        end
        i += 1
    end
    abspath(out), seed, images
end

"""
Refuse to write into a directory that already holds something other than a demo
tree, so a mistyped DEMO_DIR can't scatter synthetic files into real data.
"""
function check_target(out)
    isdir(out) || return
    entries = filter(!startswith("."), readdir(out))
    isempty(entries) || isfile(joinpath(out, "TRUTH.csv")) ||
        error("$out is not empty and has no TRUTH.csv, so it is not a demo tree; refusing to write into it")
end

function main(args = ARGS)
    out, seed, images = parse_args(args)
    check_target(out)
    # Separate streams: traces/positions must not depend on whether frames are
    # drawn, so `--no-images` yields the same .dat files `make demo` serves.
    rng = Xoshiro(seed)
    rng_img = Xoshiro(seed + 1)
    data = joinpath(out, "data")
    analysis = joinpath(out, "analysis")
    tot = joinpath(analysis, "automatic_analysis", "tot_files")
    # Remove files from an earlier generation so a re-run can't leave strays.
    for d in (data, tot)
        isdir(d) && foreach(f -> rm(joinpath(d, f)), readdir(d))
    end
    mkpath(data); mkpath(tot)
    write_setup_info(joinpath(analysis, "setup_info_20260425_181705.txt"))

    qm, dx, dy = images ? q_map() : (nothing, nothing, nothing)
    gaps = gap_mask()
    truth = String["sample,bases,phases,lattice_A,auto,expected_peaks_q,images,note"]
    scan = FIRST_SCAN
    t = DateTime(2026, 4, 26, 18, 0, 0)
    for (k, (label, nscans, comps, auto, note)) in enumerate(SAMPLES)
        # Holder 2 goes in after a >1 h break, so load segmentation finds 2 loads.
        k == LOAD_SPLIT + 1 && (t += Minute(75))
        slot = k <= LOAD_SPLIT ? k : k - LOAD_SPLIT
        model, peaks = model_trace(rng, label, comps)
        bases = String[]
        for _ in 1:nscans
            base = "$(label)_S$(scan)"; scan += 1
            push!(bases, base)
            I, sig = measure(rng, model)
            open(joinpath(tot, "$(base)_tot.dat"), "w") do io
                for j in eachindex(Q)
                    @printf(io, "%.6e %.6e %.6e\n", Q[j], I[j], sig[j])
                end
            end
            stem = base * "_0_001"
            t += Second(rand(rng, 100:160))
            # Each sample sits in its own capillary slot: HimalayaUI groups exposures
            # into samples by horizontal-position jumps (grouping.jl _cluster_slots).
            horiz = round(SLOT0_MM + SLOT_PITCH_MM * (slot - 1) + 0.05 * randn(rng); digits = 3)
            write_prp(joinpath(data, stem * ".prp"), Dates.format(t, "dd u yyyy HH:MM:SS"); horiz)
            tif = joinpath(data, stem * ".tif")
            images ? write_tif(tif, make_frame(rng_img, I, qm, dx, dy, gaps, 1.0)) : touch(tif)
        end
        phases = isempty(comps) ? "none" : join((phase_name(c[1]) for c in comps), "+")
        lattice = join((@sprintf("%g", c[2]) for c in comps), "+")
        qs = join(("$p:$q" for (p, q) in peaks), " ")
        push!(truth, join([label, join(bases, " "), phases, lattice, auto, qs,
                           string(images), csvq(note)], ","))
    end
    write(joinpath(out, "TRUTH.csv"), join(truth, "\n") * "\n")
    println("Wrote synthetic demo experiment to $out ($(length(SAMPLES)) samples, seed $seed",
            images ? "" : ", no images", ")")
    out
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
