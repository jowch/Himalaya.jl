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
#       [out_dir] [--seed 7] [--no-images] [--large [N]]
#
# `--large [N]` (default N = 300) writes a beamtime-sized experiment instead of the
# 13-sample screen: N samples cycled from SAMPLES, in holders of up to 16 slots,
# with frames per sample drawn to match real SSRL data (median 4, a tenth at 10,
# a few at 15 and 30). Frames that share a sample template and bad-frame kind are
# hard links to one TIFF, so the tree stays small on disk while every exposure
# still has its own file for the server to render and cache.
#
# Layout written (what routes_fs.jl `resolve_experiment_layout` auto-detects);
# <base> = <label>_S<scan>, one scan per exposure, 1-3 scans per sample, a few of
# them bad on purpose (BAD_SCANS):
#   <out>/data/<base>_0_001.tif        Pilatus 1M frame, Int32 counts (-1 in module gaps)
#   <out>/data/<base>_0_001.prp        SSRL PRP sidecar (energy, pipe length, position, time, beam monitors)
#   <out>/analysis/setup_info_<ts>.txt beam center + AgBe-calibrated distance
#   <out>/analysis/automatic_analysis/tot_files/<base>_tot.dat   q  I  sigma
#   <out>/TRUTH.csv                    per sample: generated phases, expected auto call,
#                                      phase-tagged reflection q values, scans to drop,
#                                      and which of those are flares (sample still there)
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
const PILATUS_MAX_COUNTS = 1_000_000
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
    ("leaked_capillary", 2, Component[],                                                 "none",           "capillary leaked before measurement; both scans are air, so the sample has no usable frame"),
]

# physics: a beamtime often runs some samples at a longer camera length to reach lower
# physics: q, so their traces stop at a lower q_max than the rest. Form-factor samples
# physics: are the usual candidates, and with no Bragg peaks none of TRUTH's q go missing.
const Q_MAX = Dict("vesicles_ULV" => 0.2)

# Bad scans a triager should drop, by sample and scan number (1 = oldest); TRUTH's
# `drop` column lists them.
const BAD_SCANS = Dict(
    "MO_H2O"           => Dict(3 => :beam_off),
    "MO_H2O_37C"       => Dict(2 => :flare),
    "MO_F127"          => Dict(1 => :flare),
    "DOPE_H2O"         => Dict(1 => :missed),
    "DOPC_H2O"         => Dict(2 => :flare),
    "MO_AMP_10"        => Dict(3 => :flare),
    "leaked_capillary" => Dict(1 => :missed, 2 => :missed),
)

# physics: beam off (shutter closed) leaves only the detector's dark floor, about
# physics: 0.01 counts per pixel in 15 s, with no q dependence.
# physics: a missed capillary leaves only air and slit (parasitic) scatter: a steep
# physics: ~q^-4 rise into the beamstop over a faint air floor, about an eighth of the
# physics: incoherent floor a filled capillary gives, so the frame is dark past the
# physics: first few hundred pixels.
const DARK_FLOOR = 0.2
const AIR_FLOOR = 4.0
bad_trace(kind::Symbol, model) =
    kind === :beam_off ? fill(DARK_FLOOR, length(Q)) :
    kind === :missed   ? @.(0.002 * Q^-4 + AIR_FLOOR) :
    kind === :flare    ? @.(FLARE_SAMPLE_FRACTION * model + flare_average(Q)) :
    error("unknown bad scan kind $kind")

# physics: the commonest bad frame is a flare: the beam clips the capillary's glass
# physics: wall. The wall scatters into a bright, narrow streak through the beamstop,
# physics: perpendicular to the (vertical) capillary axis, so horizontal on the
# physics: detector. A flat interface scatters into a rod along its normal whose
# physics: intensity falls as q^-2, so the streak reaches well past the rings.
# physics: Less of the beam crosses the sample, so its peaks weaken but stay. An
# physics: azimuthal average dilutes the streak by its share of each ring: two
# physics: Gaussian arms of half-width w on a ring of radius r cover w*sqrt(pi)/(pi r).
const FLARE_SAMPLE_FRACTION = 0.6
const FLARE_HALF_WIDTH_PX = 12.0
const FLARE_AMP = 160.0
flare_streak(q) = FLARE_AMP * q^-2
radius_px(q) = MEAN_DISTANCE_MM * 1e-3 * tan(2 * asin(q * WAVELENGTH_A / (4π))) / (PIXEL_UM * 1e-6)
flare_average(q) = flare_streak(q) * FLARE_HALF_WIDTH_PX / (sqrt(π) * radius_px(q))
const LOAD_SPLIT = 6        # samples 1..6 in holder 1, the rest in holder 2
const FIRST_SCAN = 2401
const LARGE_DEFAULT = 300
const LARGE_SLOTS = 16
const LARGE_FRAMES = [(1, 2), (2, 10), (3, 20), (4, 25), (5, 12), (6, 8), (8, 8), (10, 10), (15, 3), (30, 2)]

function draw_frames(rng)
    r = rand(rng, 1:sum(last, LARGE_FRAMES))
    for (n, w) in LARGE_FRAMES
        r <= w && return n
        r -= w
    end
    first(last(LARGE_FRAMES))
end

function sample_plan(large::Int, seed::Int)
    large == 0 && return [(label = l, template = l, nscans = n, comps = c, auto = a, note = no,
                           load = k <= LOAD_SPLIT ? 1 : 2, slot = k <= LOAD_SPLIT ? k : k - LOAD_SPLIT)
                          for (k, (l, n, c, a, no)) in enumerate(SAMPLES)]
    rng = Xoshiro(seed + 3)
    map(1:large) do k
        l, _, c, a, no = SAMPLES[mod1(k, length(SAMPLES))]
        (label = @sprintf("%s_%03d", l, k), template = l, nscans = draw_frames(rng), comps = c, auto = a,
         note = no, load = cld(k, LARGE_SLOTS), slot = mod1(k, LARGE_SLOTS))
    end
end

phase_name(::Type{P}) where {P<:Phase} = string(nameof(P))

"Linear interpolation on the global Q grid, clamped at the ends."
function interpq(x, ys)
    x <= Q[1] && return ys[1]
    x >= Q[end] && return ys[end]
    k = searchsortedlast(Q, x)
    t = (x - Q[k]) / (Q[k+1] - Q[k])
    ys[k] + t * (ys[k+1] - ys[k])
end

# physics: the flat floor is solvent and capillary scatter. Every capillary holds
# physics: the same buffer, so per unit transmitted flux the floor is the same in
# physics: every sample and buffer_blank is exactly what a subtraction should remove.
const SOLVENT_FLOOR = 30.0

"Power-law low-q upturn + flat incoherent floor, shaped like example_tot.dat."
function background(rng, formfactor::Bool)
    amp = 0.7 + 0.6 * rand(rng)
    if formfactor
        # broad bilayer form-factor bump around 0.1-0.2 Å⁻¹ on a q^-2 upturn
        return @. amp * (3.0 * Q^-2.0 + 120 * exp(-((Q - 0.14) / 0.045)^2)) + SOLVENT_FLOOR
    end
    @. amp * (0.06 * Q^-3.2 + 1.2 * Q^-1.2) + SOLVENT_FLOOR
end

# physics: traces are written as raw detector counts, as the beamline's integration
# physics: writes them, so each scales with the flux that crossed the sample during
# physics: the exposure. That flux is the incident beam (I_0, which decays between
# physics: top-ups) times the sample's transmission, and the beamstop APD downstream
# physics: records it as I_signal. Dividing a trace by its I_signal puts every
# physics: exposure on one scale; a high-q match cannot, when the sample scatters there.
# physics: The constants follow the real SSRL 2026 .prp: I_0 ≈ 1.50 V and
# physics: I_signal ≈ 7365 counts for 15 s on a lipid sample, with the APD's dark rate
# physics: (18320 /s) included in raw I_1 and removed from I_signal.
const I0_VOLTS     = 1.5007
const APD_PER_VOLT = 7365.0 / (1.5007 * 0.86)
const APD_DARK_RATE = 18320.0
# physics: I_signal is a small difference of two large numbers: raw I_1 minus the
# physics: nominal dark rate × counting time, about 97% of I_1 in the real file. The
# physics: dark rate drifts between its calibration and each exposure, so a 0.1% drift
# physics: moves I_signal by about 4%. The demo draws that drift, so the monitor scale
# physics: it shows carries the error a real one would, rather than agreeing exactly.
const APD_DARK_DRIFT = 0.001
const TRANSMISSION = Dict("buffer_blank" => 0.93, "vesicles_ULV" => 0.90)
const AIR_TRANSMISSION = 1.0

"Incident flux (relative) and transmission for one exposure."
function beam(rng, label, kind, minutes)
    i0 = (1.0 - 0.18 * mod(minutes, 120) / 120) * (1 + 0.01 * randn(rng))
    kind === :beam_off && return (i0 = 0.0, t = 0.0)
    t = kind === :missed ? AIR_TRANSMISSION : get(TRANSMISSION, label, 0.80 + 0.1 * rand(rng))
    (i0 = i0, t = t)
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

function make_frame(rng, ring, qm, dx, dy, gaps, frame_scale; flare = false)
    counts = Matrix{Int32}(undef, DET_H, DET_W)
    for i in eachindex(qm)
        # radial profile -> 2D, scaled to plausible per-pixel counts for a 15 s exposure
        λ = interpq(qm[i], ring)
        flare && (λ += flare_streak(qm[i]) * exp(-(dy[i] / FLARE_HALF_WIDTH_PX)^2))
        counts[i] = min(poisson(rng, λ * 0.05 * frame_scale), PILATUS_MAX_COUNTS)
        # beamstop + arm shadow (arm runs down-left from the beam centre)
        r = hypot(dx[i], dy[i])
        if r < 18 || (abs(dy[i] - 0.55 * dx[i]) < 6 && dx[i] < 0)
            counts[i] = 0
        end
    end
    for k in rand(rng, 1:length(counts), 6)        # a few hot pixels
        counts[k] = PILATUS_MAX_COUNTS
    end
    counts[gaps] .= -1
    counts
end

# Pilatus frames are signed 32-bit integer TIFFs. TiffImages maps that pixel type
# to Gray{Q0f31}, and image.jl reinterprets back to Int32 counts, so writing the
# raw Int32 bits as Q0f31 produces the same on-disk format as the beamline.
write_tif(path, counts::Matrix{Int32}) = save(path, Gray.(reinterpret.(Q0f31, counts)))

function write_prp(path, ts; exposure_time = 15, horiz = 58.9, i0 = 1.0, transmission = 0.86, rng = Xoshiro(0))
    stem = splitext(basename(path))[1]
    counting = exposure_time + 0.8
    transmitted = APD_PER_VOLT * I0_VOLTS * i0 * transmission * counting / 15.8 * (1 + 0.003 * randn(rng))
    dark = APD_DARK_RATE * (1 + APD_DARK_DRIFT * randn(rng))
    i_1 = round(transmitted + dark * counting; digits = 1)
    i_signal = round(i_1 - APD_DARK_RATE * counting; digits = 1)
    open(path, "w") do io
        println(io, "Image file name: $stem.tif")
        println(io, "Time this file was written: $ts")
        println(io, "I_Signal (i.e. I1 in header file) is from: i_apd")
        println(io, "I_signal=$i_signal (corrected for dark counts)")
        @printf(io, "I_0(daqIon3.0)=%.4f\n", I0_VOLTS * i0)
        println(io, "I_1(daqCounter1.0)=$i_1")
        println(io, "I_apd(daqCounter1.0)=$i_1")
        println(io, "Exposure time=$exposure_time")
        println(io, "Counting time=$counting")
        println(io, "i1 dark counts=$APD_DARK_RATE per sec")
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
    seed = 7; images = true; large = 0
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--seed"
            i < length(args) || error("--seed needs a value")
            seed = parse(Int, args[i+1]); i += 1
        elseif a == "--no-images"
            images = false
        elseif a == "--large"
            if i < length(args) && !startswith(args[i+1], "-") && tryparse(Int, args[i+1]) !== nothing
                large = parse(Int, args[i+1]); i += 1
            else
                large = LARGE_DEFAULT
            end
            large > 0 || error("--large needs a positive sample count")
        elseif startswith(a, "-")
            error("unknown option $a (usage: make_demo_data.jl [out_dir] [--seed N] [--no-images] [--large [N]])")
        else
            out = a
        end
        i += 1
    end
    abspath(out), seed, images, large
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
    out, seed, images, large = parse_args(args)
    plan = sample_plan(large, seed)
    shared_tif = Dict{Tuple{String,Union{Symbol,Nothing}},String}()
    check_target(out)
    # Separate streams: traces/positions must not depend on whether frames are
    # drawn, so `--no-images` yields the same .dat files `make demo` serves.
    rng = Xoshiro(seed)
    rng_img = Xoshiro(seed + 1)
    rng_beam = Xoshiro(seed + 2)
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
    truth = String["sample,bases,phases,lattice_A,auto,expected_peaks_q,images,drop,flare,note"]
    scan = FIRST_SCAN
    t = DateTime(2026, 4, 26, 18, 0, 0)
    for (k, smp) in enumerate(plan)
        (; label, template, nscans, comps, auto, note, slot) = smp
        # Each holder after the first goes in after a >1 h break, so load segmentation
        # finds one load per holder.
        k > 1 && smp.load != plan[k-1].load && (t += Minute(75))
        model, peaks = model_trace(rng, template, comps)
        bad = get(BAD_SCANS, template, Dict{Int,Symbol}())
        bases = String[]; drops = String[]; flares = String[]
        for n in 1:nscans
            base = "$(label)_S$(scan)"; scan += 1
            push!(bases, base)
            kind = get(bad, n, nothing)
            kind === nothing || push!(drops, base)
            kind === :flare && push!(flares, base)
            I, sig = measure(rng, kind === nothing ? model : bad_trace(kind, model))
            b = beam(rng_beam, template, kind, Dates.value(t - DateTime(2026, 4, 26, 18, 0, 0)) / 60_000)
            flux = kind === :beam_off ? 1.0 : b.i0 * b.t / 0.86
            ring = kind === :flare ? FLARE_SAMPLE_FRACTION .* model : I
            I = I .* flux; sig = sig .* flux
            open(joinpath(tot, "$(base)_tot.dat"), "w") do io
                for j in eachindex(Q)
                    Q[j] <= get(Q_MAX, template, Inf) || continue
                    @printf(io, "%.6e %.6e %.6e\n", Q[j], I[j], sig[j])
                end
            end
            stem = base * "_0_001"
            t += Second(rand(rng, 100:160))
            # Each sample sits in its own capillary slot: HimalayaUI groups exposures
            # into samples by horizontal-position jumps (grouping.jl _cluster_slots).
            horiz = round(SLOT0_MM + SLOT_PITCH_MM * (slot - 1) + 0.05 * randn(rng); digits = 3)
            write_prp(joinpath(data, stem * ".prp"), Dates.format(t, "dd u yyyy HH:MM:SS");
                      horiz, i0 = b.i0, transmission = b.t, rng = rng_beam)
            tif = joinpath(data, stem * ".tif")
            if !images
                touch(tif)
            elseif large > 0 && haskey(shared_tif, (template, kind))
                hardlink(shared_tif[(template, kind)], tif)
            else
                write_tif(tif, make_frame(rng_img, ring, qm, dx, dy, gaps, flux; flare = kind === :flare))
                shared_tif[(template, kind)] = tif
            end
        end
        phases = isempty(comps) ? "none" : join((phase_name(c[1]) for c in comps), "+")
        lattice = join((@sprintf("%g", c[2]) for c in comps), "+")
        qs = join(("$p:$q" for (p, q) in peaks), " ")
        push!(truth, join([label, join(bases, " "), phases, lattice, auto, qs,
                           string(images), join(drops, " "), join(flares, " "), csvq(note)], ","))
    end
    write(joinpath(out, "TRUTH.csv"), join(truth, "\n") * "\n")
    println("Wrote synthetic demo experiment to $out ($(length(plan)) samples, seed $seed",
            images ? "" : ", no images", ")")
    out
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
