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
# `--large [--experiments 10] [--samples 300]` instead writes <out>/large_NN/ trees at
# realistic scale (`make demo-large-data`); see docs/demo.md "Large mode". There each
# .tif is a hardlink to one of ~24 pooled frames, not rendered from its own trace.
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

using Himalaya: Phase, Lamellar, Hexagonal, Square, Pn3m, Im3m, Ia3d, phaseratios
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
    Square    => [1.0, 0.40, 0.20, 0.12, 0.06, 0.04, 0.03, 0.02, 0.015, 0.01],
    Pn3m      => [1.0, 0.75, 0.12, 0.10, 0.06, 0.05, 0.03, 0.02, 0.015],
    Im3m      => [1.0, 0.55, 0.30, 0.10, 0.06, 0.035, 0.02],
    Ia3d      => [1.0, 0.65, 0.06, 0.08, 0.04, 0.03, 0.02],
)

"q of the first allowed reflection for lattice parameter `a` (Å)."
first_peak_q(::Type{Lamellar}, a) = 2π / a                    # a = d-spacing
first_peak_q(::Type{Hexagonal}, a) = 4π / (sqrt(3) * a)       # a = rod spacing
first_peak_q(::Type{Square}, a) = 2π / a                      # a = cell edge; (1,0) reflection
first_peak_q(::Type{P}, a) where {P<:Phase} = 2π * first(phaseratios(P)) / a  # cubic: 2π√N/a

# (phase, lattice parameter Å, peak height relative to background at q1, peak sigma)
# Heights sit in the real range (on SSRL lipid traces the strongest indexed peak
# is a median ~0.4× its local background, ~8× at the 90th percentile) and widths
# near it (half-maximum width ~0.003-0.010 Å⁻¹). The two-phase samples use the
# sharp end, since broader peaks merge their neighbouring first reflections.
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
    ("MO_H2O",         3, Component[(Pn3m, 102.0, 5.0, 2.1e-3)],                          "Pn3m",           "monoolein in excess water, clean Pn3m"),
    ("MO_H2O_37C",     2, Component[(Pn3m, 95.0, 4.0, 2.1e-3)],                           "Pn3m",           "Pn3m, lattice shrinks on heating"),
    ("MO_PEG5",        2, Component[(Pn3m, 108.0, 3.0, 1.6e-3), (Im3m, 138.0, 2.0, 1.6e-3)], "Im3m",          "Pn3m + Im3m coexistence (Bonnet ratio ~1.28). Pn3m √6 overlaps Im3m √10, so auto_group keeps Im3m; Pn3m is the second hypothesis to promote by hand"),
    ("MO_F127",        2, Component[(Im3m, 145.0, 3.0, 2.4e-3)],                          "Im3m",           "polymer-stabilised Im3m"),
    ("MO_lowH2O",      2, Component[(Ia3d, 150.0, 4.0, 2.1e-3)],                          "Ia3d",           "low hydration Ia3d"),
    ("DOPE_H2O",       3, Component[(Hexagonal, 74.0, 6.0, 2.7e-3)],                      "Hexagonal",      "DOPE inverse hexagonal"),
    ("DOPC_H2O",       2, Component[(Lamellar, 63.0, 8.0, 3.0e-3)],                     "Lamellar",       "DOPC multilamellar"),
    ("DOPC_DOPE_1to1", 2, Component[(Lamellar, 60.0, 4.0, 1.3e-3), (Hexagonal, 72.0, 4.0, 1.3e-3)], "Hexagonal+Lamellar", "La + HII coexistence"),
    ("MO_AMP_weak",    2, Component[(Pn3m, 99.0, 0.12, 2.7e-3)],                            "none",           "weak Pn3m near the noise floor; too weak to index automatically"),
    ("MO_AMP_10",      3, Component[(Pn3m, 96.0, 6.0, 1.0e-3), (Hexagonal, 53.0, 6.0, 1.0e-3)], "Hexagonal+Pn3m", "Pn3m drifting to HII under peptide"),
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
# physics: a missed capillary leaves only air and slit (parasitic) scatter: half the
# physics: low-q rise a filled capillary shows, over a faint air floor about an eighth
# physics: of the solvent floor, so the frame is dim past the first hundred pixels.
const DARK_FLOOR = 0.002
bad_trace(kind::Symbol, model) =
    kind === :beam_off ? fill(DARK_FLOOR, length(Q)) :
    kind === :missed   ? @.(AIR_FLOOR + 0.5 * SOLVENT_FLOOR * RISE * Q^-RISE_EXP) :
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
const FLARE_AMP = 90.0
flare_streak(q) = FLARE_AMP * q^-2
radius_px(q) = MEAN_DISTANCE_MM * 1e-3 * tan(2 * asin(q * WAVELENGTH_A / (4π))) / (PIXEL_UM * 1e-6)
flare_average(q) = flare_streak(q) * FLARE_HALF_WIDTH_PX / (sqrt(π) * radius_px(q))
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

# physics: the flat floor is solvent and capillary scatter. Every capillary holds
# physics: the same buffer, so per unit transmitted flux the floor is the same in
# physics: every sample and buffer_blank is exactly what a subtraction should remove.
# physics: Over it sit a steep rise into the beamstop (parasitic and large-scale
# physics: scatter) and a broad bilayer form-factor hump near 0.11 Å⁻¹. The shape
# physics: follows the median of SSRL 1M lipid traces: about 100× the floor at
# physics: q = 0.006, 2× at 0.03, a hump of about 2× at 0.12, flat past 0.25.
const SOLVENT_FLOOR = 17.0
const RISE, RISE_EXP = 4.5e-5, 2.8
const HUMP, HUMP_Q, HUMP_W = 1.1, 0.11, 0.045
# physics: a trace is the azimuthal mean per pixel in units where the solvent floor
# physics: is ~100 counts per pixel on a 15 s frame, as on the beamline.
const PX_PER_TRACE = 6.0
const AIR_FLOOR = SOLVENT_FLOOR / 8

"Solvent floor × (1 + low-q rise + bilayer hump), with the rise and hump varying by sample."
function background(rng, formfactor::Bool)
    rise = RISE * exp(0.5 * randn(rng))
    hump = (formfactor ? 3.0 : HUMP) * (0.6 + 0.8 * rand(rng))
    @. SOLVENT_FLOOR * (1 + rise * Q^-RISE_EXP + hump * exp(-0.5 * ((Q - HUMP_Q) / HUMP_W)^2))
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

# physics: each trace point averages every detector pixel in its q bin, so its
# physics: counting error is sqrt(counts per pixel × pixels in the bin) over that
# physics: total: about 0.3% at high q, as on the beamline. NPIX counts the live
# physics: pixels per Q bin from the detector geometry (gaps included, which is close
# physics: enough for noise).
const NPIX = let
    n = zeros(length(Q)); dq = Q[2] - Q[1]
    for c in 1:DET_W, r in 1:DET_H
        rm = hypot(c - 1 - BEAM_X, DET_H - r - BEAM_Y) * PIXEL_UM * 1e-6
        q = 4π * sin(atan(rm / (MEAN_DISTANCE_MM * 1e-3)) / 2) / WAVELENGTH_A
        k = round(Int, (q - Q[1]) / dq) + 1
        1 <= k <= length(Q) && (n[k] += 1)
    end
    max.(n, 1.0)
end

"One noisy measurement of a model trace; σ is the noise actually applied."
function measure(rng, I)
    sig = @. sqrt(max(I, 1e-6) / (PX_PER_TRACE * NPIX))
    @. max(I + randn(rng) * sig, 1e-6), sig
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
function q_map(distance_mm = MEAN_DISTANCE_MM)
    qm = Matrix{Float64}(undef, DET_H, DET_W)
    dx = similar(qm); dy = similar(qm)
    for c in 1:DET_W, r in 1:DET_H
        x = c - 1; yup = DET_H - r
        dx[r, c] = x - BEAM_X; dy[r, c] = yup - BEAM_Y
        rm = hypot(dx[r, c], dy[r, c]) * PIXEL_UM * 1e-6
        tth = atan(rm / (distance_mm * 1e-3))
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

const BEAMSTOP_PX = 10.0

function make_frame(rng, ring, qm, dx, dy, gaps, frame_scale; flare = false)
    counts = Matrix{Int32}(undef, DET_H, DET_W)
    for i in eachindex(qm)
        # radial profile -> 2D, scaled to plausible per-pixel counts for a 15 s exposure
        λ = interpq(qm[i], ring)
        flare && (λ += flare_streak(qm[i]) * exp(-(dy[i] / FLARE_HALF_WIDTH_PX)^2))
        counts[i] = min(poisson(rng, λ * PX_PER_TRACE * frame_scale), PILATUS_MAX_COUNTS)
        # beamstop + arm shadow (arm runs down-left from the beam centre)
        r = hypot(dx[i], dy[i])
        if r < BEAMSTOP_PX || (abs(dy[i] - 0.55 * dx[i]) < 6 && dx[i] < 0)
            counts[i] = 0
        end
    end
    counts[gaps] .= -1
    counts
end

# Pilatus frames are signed 32-bit integer TIFFs. TiffImages maps that pixel type
# to Gray{Q0f31}, and image.jl reinterprets back to Int32 counts, so writing the
# raw Int32 bits as Q0f31 produces the same on-disk format as the beamline.
write_tif(path, counts::Matrix{Int32}) = save(path, Gray.(reinterpret.(Q0f31, counts)))

# physics: with the shutter closed the ion chamber's electrometer reads its offset, a
# physics: few mV of either sign (the idle Vt channel on the real 2026 .prp reads
# physics: -0.0045 V), not 0. Drawn from the file stem so the shared beam stream
# physics: and every other frame's output stay as they were.
function electrometer_offset(stem)
    r = Xoshiro(hash(stem))
    (rand(r, Bool) ? 1 : -1) * (0.001 + 0.004 * rand(r))
end

function write_prp(path, ts; exposure_time = 15, horiz = 58.9, pipe_mm = PIPE_LENGTH_MM, i0 = 1.0, transmission = 0.86, rng = Xoshiro(0))
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
        @printf(io, "I_0(daqIon3.0)=%.4f\n", i0 == 0 ? electrometer_offset(stem) : I0_VOLTS * i0)
        println(io, "I_1(daqCounter1.0)=$i_1")
        println(io, "I_apd(daqCounter1.0)=$i_1")
        println(io, "Exposure time=$exposure_time")
        println(io, "Counting time=$counting")
        println(io, "i1 dark counts=$APD_DARK_RATE per sec")
        println(io, "Beam energy=$ENERGY_EV eV")
        println(io, "Pipe length=$pipe_mm mm")
        println(io, "Scan motor=motor\nScan range=delta mm")
        println(io, "Horizontal position=$horiz mm")
        println(io, "Vertical position=26.55 mm\ndispx position=9.8 mm\ndispy position=35.9 mm")
        println(io, "Detector=Pilatus 1M")
        println(io, "Detector mode=Pilatus; 0 for normal and 1 for dezingered")
        println(io, "Phi position=157\nAttenuator=0 um Al")
    end
end

function write_setup_info(path; distance_mm = MEAN_DISTANCE_MM)
    open(path, "w") do io
        print(io, "File used for calibration:\n../data/agbe_1p7m_S1958_0_001.tif\n\n")
        print(io, "Beam center is at ($BEAM_X, $BEAM_Y)\n\n")
        println(io, "Calculated detector distances are:")
        for (i, d) in enumerate([1806.0, 1811.1, 1811.5])
            @printf(io, "Ring %d:\t pixel %.2f:\t %s mm\n", i, 247.82 * i, d)
        end
        println(io, "-"^58)
        println(io, "Mean distance:\t\t $distance_mm mm")
    end
end

csvq(s) = occursin(r"[,\"]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s

const USAGE = "usage: make_demo_data.jl [out_dir] [--seed N] [--no-images] " *
              "[--large [--experiments N] [--samples M]]"

function parse_args(args)
    out = joinpath(dirname(dirname(dirname(@__DIR__))), "scratch", "demo-experiment")
    o = (seed = 7, images = true, large = false, experiments = 10, samples = 300)
    i = 1
    intopt(name) = (i < length(args) || error("$name needs a value"); parse(Int, args[i+1]))
    while i <= length(args)
        a = args[i]
        if a == "--seed"
            o = merge(o, (seed = intopt(a),)); i += 1
        elseif a == "--experiments"
            o = merge(o, (experiments = intopt(a),)); i += 1
        elseif a == "--samples"
            o = merge(o, (samples = intopt(a),)); i += 1
        elseif a == "--no-images"
            o = merge(o, (images = false,))
        elseif a == "--large"
            o = merge(o, (large = true,))
        elseif startswith(a, "-")
            error("unknown option $a ($USAGE)")
        else
            out = a
        end
        i += 1
    end
    merge(o, (out = abspath(out),))
end

"""
Refuse to write into a directory that already holds something other than a demo
tree, so a mistyped DEMO_DIR can't scatter synthetic files into real data. A large
tree is a directory of experiment directories, each with its own TRUTH.csv.
"""
const LARGE_MARKER = ".demo-large-complete"
is_large_entry(e) = occursin(r"^large_\d+$", e)

function check_target(out, large::Bool = false)
    isdir(out) || return
    entries = filter(!startswith("."), readdir(out))
    isempty(entries) && return
    curated = isfile(joinpath(out, "TRUTH.csv"))
    largetree = all(is_large_entry, entries)
    curated && large && error("$out holds the curated demo; write --large somewhere else")
    largetree && !large && error("$out holds a --large demo; write the curated demo somewhere else")
    curated || largetree ||
        error("$out is not empty and is not a demo tree; refusing to write into it")
end

const TRUTH_HEADER = "sample,bases,phases,lattice_A,auto,expected_peaks_q,images,drop,flare,note"
const T0_DEFAULT = DateTime(2026, 4, 26, 18, 0, 0)

"""
Reset an experiment tree and return its (data, analysis, tot) directories. Files
from an earlier generation are removed so a re-run can't leave strays.
"""
function prepare_tree(out)
    data = joinpath(out, "data")
    analysis = joinpath(out, "analysis")
    tot = joinpath(analysis, "automatic_analysis", "tot_files")
    for d in (data, tot)
        isdir(d) && foreach(f -> rm(joinpath(d, f)), readdir(d))
    end
    mkpath(data); mkpath(tot)
    data, analysis, tot
end

"""
Write one exposure's `.dat`, `.prp` and `.tif` and return the time after it.
`emit_tif(path, ring, flux)` writes the frame; `ring` is the radial profile a
rendered frame would show. RNG draws happen in a fixed order (measure, beam,
exposure gap, slot jitter, monitors) that the default demo's output depends on.
"""
function write_scan!(rng, rng_beam, dirs, base, label, kind, model, t, t0, slot, emit_tif;
                     pipe_mm = PIPE_LENGTH_MM)
    I, sig = measure(rng, kind === nothing ? model : bad_trace(kind, model))
    b = beam(rng_beam, label, kind, Dates.value(t - t0) / 60_000)
    flux = kind === :beam_off ? 1.0 : b.i0 * b.t / 0.86
    ring = kind === :flare ? FLARE_SAMPLE_FRACTION .* model : I
    I = I .* flux; sig = sig .* flux
    open(joinpath(dirs.tot, "$(base)_tot.dat"), "w") do io
        for j in eachindex(Q)
            Q[j] <= get(Q_MAX, label, Inf) || continue
            @printf(io, "%.6e %.6e %.6e\n", Q[j], I[j], sig[j])
        end
    end
    stem = base * "_0_001"
    t += Second(rand(rng, 100:160))
    # Each sample sits in its own capillary slot: HimalayaUI groups exposures
    # into samples by horizontal-position jumps (grouping.jl _cluster_slots).
    horiz = round(SLOT0_MM + SLOT_PITCH_MM * (slot - 1) + 0.05 * randn(rng); digits = 3)
    write_prp(joinpath(dirs.data, stem * ".prp"), Dates.format(t, "dd u yyyy HH:MM:SS");
              horiz, pipe_mm, i0 = b.i0, transmission = b.t, rng = rng_beam)
    emit_tif(joinpath(dirs.data, stem * ".tif"), ring, flux, kind === :flare)
    t
end

truth_row(label, bases, comps, peaks, auto, images, drops, flares, note) =
    join([label, join(bases, " "),
          isempty(comps) ? "none" : join((phase_name(c[1]) for c in comps), "+"),
          join((@sprintf("%g", c[2]) for c in comps), "+"), auto,
          join(("$p:$q" for (p, q) in peaks), " "),
          string(images), join(drops, " "), join(flares, " "), csvq(note)], ",")

function main_curated(o)
    out, seed, images = o.out, o.seed, o.images
    # Separate streams: traces/positions must not depend on whether frames are
    # drawn, so `--no-images` yields the same .dat files `make demo` serves.
    rng = Xoshiro(seed)
    rng_img = Xoshiro(seed + 1)
    rng_beam = Xoshiro(seed + 2)
    data, analysis, tot = prepare_tree(out)
    dirs = (data = data, tot = tot)
    write_setup_info(joinpath(analysis, "setup_info_20260425_181705.txt"))

    qm, dx, dy = images ? q_map() : (nothing, nothing, nothing)
    gaps = gap_mask()
    emit_tif(tif, ring, flux, flare) =
        images ? write_tif(tif, make_frame(rng_img, ring, qm, dx, dy, gaps, flux; flare)) : touch(tif)
    truth = String[TRUTH_HEADER]
    scan = FIRST_SCAN
    t = T0_DEFAULT
    for (k, (label, nscans, comps, auto, note)) in enumerate(SAMPLES)
        # Holder 2 goes in after a >1 h break, so load segmentation finds 2 loads.
        k == LOAD_SPLIT + 1 && (t += Minute(75))
        slot = k <= LOAD_SPLIT ? k : k - LOAD_SPLIT
        model, peaks = model_trace(rng, label, comps)
        bad = get(BAD_SCANS, label, Dict{Int,Symbol}())
        bases = String[]; drops = String[]; flares = String[]
        for n in 1:nscans
            base = "$(label)_S$(scan)"; scan += 1
            push!(bases, base)
            kind = get(bad, n, nothing)
            kind === nothing || push!(drops, base)
            kind === :flare && push!(flares, base)
            t = write_scan!(rng, rng_beam, dirs, base, label, kind, model, t, T0_DEFAULT, slot, emit_tif)
        end
        push!(truth, truth_row(label, bases, comps, peaks, auto, images, drops, flares, note))
    end
    write(joinpath(out, "TRUTH.csv"), join(truth, "\n") * "\n")
    println("Wrote synthetic demo experiment to $out ($(length(SAMPLES)) samples, seed $seed",
            images ? "" : ", no images", ")")
    out
end

# ---- large mode ----------------------------------------------------------------
# Realistic scale for performance work: several experiments of hundreds of samples.
# Mix targets are aggregate figures from real beamtimes (no real names or values).

# Phase weights among samples that have peaks (renormalised on use).
const LARGE_PHASE_WEIGHTS = [(Lamellar, 0.27), (Im3m, 0.16), (Square, 0.13), (Pn3m, 0.06),
                             (Hexagonal, 0.03), (Ia3d, 0.015)]
# Typical lattice parameter (Å) per phase; each sample spreads ±15% around it.
const LARGE_LATTICE = Dict{Type{<:Phase},Float64}(Lamellar => 63.0, Hexagonal => 74.0, Square => 80.0,
                                                  Pn3m => 100.0, Im3m => 140.0, Ia3d => 150.0)
const LARGE_FLIGHT_PATHS_MM = (800, 1700)
const NO_PEAK_FRACTION = 0.40
const TWO_PHASE_FRACTION = 0.15
const BAD_FRAME_FRACTION = 0.06
const POOL_SIZE = 24

"Draw from `pairs` of (value, weight)."
function weighted(rng, pairs)
    r = rand(rng) * sum(last, pairs)
    for (v, w) in pairs
        (r -= w) <= 0 && return v
    end
    first(last(pairs))
end

"Frames per sample: 4 mostly, with a tail of short, long and very long series."
function draw_nframes(rng)
    r = rand(rng)
    r < 0.82 && return 4
    r < 0.91 && return 2
    r < 0.93 && return 1
    r < 0.96 && return 5
    r < 0.99 && return rand(rng, 8:10)
    rand(rng, 15:30)
end

# Heights: median 0.4, 90th percentile 8 (σ_ln = ln(20)/1.2816). Sigmas: median
# 3.4e-3, 10th percentile 1.2e-3 (σ_ln = ln(3.4/1.2)/1.2816).
# Clamped so the log-normal tails stay physical: real strongest peaks reach ~8× their
# background at the 90th percentile, and half-max widths stay under ~0.03 Å⁻¹.
draw_height(rng) = clamp(0.4 * exp(2.337 * randn(rng)), 0.02, 20.0)
draw_sigma(rng) = clamp(3.4e-3 * exp(0.8127 * randn(rng)), 8e-4, 0.012)

function draw_component(rng, P)
    a = LARGE_LATTICE[P] * (1 + 0.15 * (2rand(rng) - 1))
    (P, a, draw_height(rng), draw_sigma(rng))::Component
end

"Components of one large-mode sample: none, one phase or two distinct phases."
function draw_components(rng)
    rand(rng) < NO_PEAK_FRACTION && return Component[]
    P1 = weighted(rng, LARGE_PHASE_WEIGHTS)
    comps = Component[draw_component(rng, P1)]
    if rand(rng) < TWO_PHASE_FRACTION
        P2 = weighted(rng, filter(p -> first(p) != P1, LARGE_PHASE_WEIGHTS))
        push!(comps, draw_component(rng, P2))
    end
    comps
end

"Bad kind for one frame, or `nothing`: mostly flares, some missed capillaries and beam-off."
function draw_bad(rng)
    rand(rng) < BAD_FRAME_FRACTION || return nothing
    weighted(rng, [(:flare, 0.7), (:missed, 0.2), (:beam_off, 0.1)])
end

"Slot count of the next holder."
draw_holder_size(rng) = rand(rng, 10:16)

"""
Pool of rendered frames: no-peak (some with the vesicle hump), weak and strong
examples of each phase, plus flare, missed and beam-off frames. Returns
`Dict(key => path)` with keys `:none`, `:flare`, `:missed`, `:beam_off` and phase
names, each mapping to a vector of paths.
"""
function render_pool(rng, rng_img, dir, qm, dx, dy, gaps)
    mkpath(dir)
    pool = Dict{Symbol,Vector{String}}()
    n = 0
    function add!(key, label, comps; kind = nothing)
        model, _ = model_trace(rng, label, comps)
        ring = kind === :flare ? FLARE_SAMPLE_FRACTION .* model :
               kind === nothing ? model : bad_trace(kind, model)
        counts = make_frame(rng_img, ring, qm, dx, dy, gaps, 1.0; flare = kind === :flare)
        path = joinpath(dir, @sprintf("pool_%02d.tif", n += 1))
        write_tif(path, counts)
        push!(get!(pool, key, String[]), path)
    end
    strong(P) = Component[(P, LARGE_LATTICE[P], 3.0, 2.5e-3)]
    weak(P) = Component[(P, LARGE_LATTICE[P], 0.25, 3.4e-3)]
    for P in first.(LARGE_PHASE_WEIGHTS)
        add!(Symbol(phase_name(P)), "pool", strong(P))
        add!(Symbol(phase_name(P)), "pool", weak(P))
    end
    for i in 1:5
        add!(:none, i <= 2 ? "vesicles_pool" : "pool", Component[])
    end
    add!(:flare, "pool", strong(Lamellar); kind = :flare)
    add!(:flare, "pool", Component[]; kind = :flare)
    add!(:missed, "pool", Component[]; kind = :missed)
    add!(:beam_off, "pool", Component[]; kind = :beam_off)
    for P in (Im3m, Lamellar, Square)
        add!(Symbol(phase_name(P)), "pool", Component[(P, LARGE_LATTICE[P] * 1.08, 1.0, 3.0e-3)])
    end
    @assert n == POOL_SIZE
    pool
end

"Hardlink `src` to `dst`, copying where the filesystem refuses."
function hardlink(src, dst)
    try
        Base.hardlink(src, dst)
    catch
        cp(src, dst)
    end
end

"Pool frame for an exposure: matched to its failure kind, else to the sample's first phase."
function pool_frame(rng, pool, kind, comps)
    key = kind !== nothing ? kind :
          isempty(comps) ? :none : Symbol(phase_name(first(comps)[1]))
    rand(rng, pool[key])
end

"One large-mode experiment `e` of `nsamples` samples written under `out`."
function write_large_experiment(out, e, nsamples, seed, images)
    # A stream per experiment, so experiment e is the same whatever N is.
    rng = Xoshiro(hash((seed, e)))
    rng_img = Xoshiro(hash((seed, e, :img)))
    rng_beam = Xoshiro(hash((seed, e, :beam)))
    rng_pick = Xoshiro(hash((seed, e, :pick)))
    pipe_mm = rand(rng, LARGE_FLIGHT_PATHS_MM)
    distance_mm = MEAN_DISTANCE_MM - (PIPE_LENGTH_MM - pipe_mm)
    t0 = DateTime(2026, 1, 5, 9, 0, 0) + Day(9 * (e - 1) + rand(rng, 0:3))
    tag = @sprintf("L%02d", e)
    data, analysis, tot = prepare_tree(out)
    dirs = (data = data, tot = tot)
    write_setup_info(joinpath(analysis, "setup_info_$(Dates.format(t0, "yyyymmdd"))_090000.txt"); distance_mm)

    pooldir = joinpath(out, ".pool")
    pool = nothing
    if images
        qm, dx, dy = q_map(distance_mm)
        pool = render_pool(rng, rng_img, pooldir, qm, dx, dy, gap_mask())
    end
    emit_tif(pk) = (tif, _, _, _) -> images ? hardlink(pk(), tif) : touch(tif)

    truth = String[TRUTH_HEADER]
    scan = FIRST_SCAN
    t = t0
    slot = 0; holder_size = draw_holder_size(rng)
    for k in 1:nsamples
        if slot == holder_size                       # next holder after a >1 h break
            t += Minute(rand(rng, 70:120)); slot = 0; holder_size = draw_holder_size(rng)
        end
        slot += 1
        label = @sprintf("%s_N%04d", tag, k)
        comps = draw_components(rng)
        model, peaks = model_trace(rng, label, comps)
        bases = String[]; drops = String[]; flares = String[]
        for _ in 1:draw_nframes(rng)
            base = "$(label)_S$(scan)"; scan += 1
            push!(bases, base)
            kind = draw_bad(rng)
            kind === nothing || push!(drops, base)
            kind === :flare && push!(flares, base)
            pk = () -> pool_frame(rng_pick, pool, kind, comps)
            t = write_scan!(rng, rng_beam, dirs, base, label, kind, model, t, t0, slot, emit_tif(pk); pipe_mm)
        end
        note = isempty(comps) ? "no Bragg peaks" : "generated phases"
        push!(truth, truth_row(label, bases, comps, peaks, "", images, drops, flares, note))
    end
    isdir(pooldir) && rm(pooldir; recursive = true)
    write(joinpath(out, "TRUTH.csv"), join(truth, "\n") * "\n")
    (samples = nsamples, pipe_mm)
end

function main_large(o)
    mkpath(o.out)
    isfile(joinpath(o.out, LARGE_MARKER)) && Base.Filesystem.unlink(joinpath(o.out, LARGE_MARKER))
    for e in 1:o.experiments
        name = @sprintf("large_%02d", e)
        r = write_large_experiment(joinpath(o.out, name), e, o.samples, o.seed, o.images)
        println("  $name: $(r.samples) samples, flight path $(r.pipe_mm) mm")
    end
    touch(joinpath(o.out, LARGE_MARKER))
    println("Wrote $(o.experiments) large synthetic experiments to $(o.out) (seed $(o.seed)",
            o.images ? ", images shared from a pool" : ", no images", ")")
    o.out
end

function main(args = ARGS)
    o = parse_args(args)
    check_target(o.out, o.large)
    o.large ? main_large(o) : main_curated(o)
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
