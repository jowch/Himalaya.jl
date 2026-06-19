using DelimitedFiles
using Statistics: median

# Hand-labeled peak positions in three real SAXS integration traces.
# Labels were placed via free-form clicking on data points (not snapped to
# local maxima), so 1-2 sample offsets between label and detected peak are
# normal.
const LABELED_PEAKS = Dict(
    "example_tot.dat" => [
        0.053439, 0.065558, 0.075081, 0.091962, 0.107111,
        0.113171, 0.124858, 0.130052, 0.139142, 0.155590,
    ],
    "cubic_tot.dat" => [
        0.046115, 0.058678, 0.065177, 0.071675, 0.079906,
        0.082939, 0.092037, 0.101568, 0.112832, 0.117164,
        0.121496, 0.124529, 0.130594, 0.137959, 0.143591,
        0.145324, 0.152689, 0.159620, 0.171318, 0.177816,
        0.189513, 0.194712, 0.200777,
    ],
    "form-factor_tot.dat" => Float64[],
)

# Regression floors. Numbers reflect current algorithm behavior on these
# specific traces; tightening these (raising recall, lowering spurious)
# requires algorithm work and should bump these floors with the same commit.
# See docs/peak-finding.md for
# the design and the known-hard-case discussion.
const RECALL_FLOOR = Dict(
    # Lowered 7→6 when SG edges moved to the "interp" boundary scheme. The old
    # 7 was a kneedle accident: a post-trimmed high-q noise candidate (idx 921)
    # held a near-zero sharpness under mirror padding, anchoring the bottom of
    # the sorted-sharpness curve and dragging the elbow below the natural
    # 5.4→1.6 cliff. That admitted a real but weak labeled peak (idx 274,
    # sharpness 1.61) whose recovery depended on the padding artifact rather
    # than on clearing the honest threshold. Interp gives idx 921 its honest
    # value, the elbow snaps to the cliff (5.40), and 6/10 is the principled
    # recall (idx 274 now falls below the floor; spurious stays 0).
    "example_tot.dat"     => 6,   # current: 6/10
    "cubic_tot.dat"       => 8,   # current: 8/23 — below production quality; TODO: iterate on sharpness method / CWT scale bank / normalize_by_σ default (see spec §Open Issues)
    "form-factor_tot.dat" => 0,   # no peaks expected
)
const SPURIOUS_CEILING = Dict(
    "example_tot.dat"     => 0,   # current: 0 spurious
    "cubic_tot.dat"       => 0,   # current: 0 spurious
    "form-factor_tot.dat" => 0,   # current: 0 spurious — eliminated by
                                  # high-q trim + relative-prominence floor
)

function load_trace(name)
    A = readdlm(joinpath(@__DIR__, "data", name))
    A[:, 1], A[:, 2], A[:, 3]
end

@testset "findpeaks on $name" for name in keys(LABELED_PEAKS)
    q, I, σ = load_trace(name)
    pk = findpeaks(q, I, σ)
    expected = LABELED_PEAKS[name]
    tol = 2 * median(diff(q))

    recovered = sum(any(abs.(pk.q .- q_label) .< tol) for q_label in expected; init = 0)
    spurious = sum(!any(abs.(q_found .- q_label) .< tol for q_label in expected) for q_found in pk.q; init = 0)

    @test recovered >= RECALL_FLOOR[name]
    @test spurious <= SPURIOUS_CEILING[name]
end
