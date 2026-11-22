using DelimitedFiles
using Random
using Statistics: median

# Tests for findpeaks(...; method = :hybrid). See docs/peak-finding.md.

const HYB_N = 970
const HYB_Q = collect(range(0.01, 0.4; length = HYB_N))
const HYB_DQ = step(range(0.01, 0.4; length = HYB_N))

# Power law + constant + narrow peaks (strong, weak, shoulder) + Poisson-like noise.
function _hybrid_trace(; seed = 1, peaks = true, extra = (), strong = [0.07, 0.115, 0.16, 0.23, 0.3])
    rng = MersenneTwister(seed)
    q = HYB_Q
    I = 6000.0 .* (q ./ 0.04) .^ (-4) .+ 30.0
    weak = [0.34]
    shoulder = [0.04]           # sits on the steep low-q flank
    if peaks
        for (q0, h, γ) in vcat([(x, 400.0, 0.0012) for x in strong],
                               [(x, 140.0, 0.0012) for x in weak],
                               [(x, 1100.0, 0.0012) for x in shoulder])
            I .+= h ./ (1 .+ ((q .- q0) ./ γ) .^ 2)
        end
    else
        I .+= 150.0 .* exp.(-((q .- 0.2) ./ 0.08) .^ 2)   # broad hump only
    end
    for (q0, h, γ) in extra
        I .+= h ./ (1 .+ ((q .- q0) ./ γ) .^ 2)
    end
    σ = sqrt.(I)
    I .+= σ .* randn(rng, length(q))
    q, I, σ, (strong = strong, weak = weak, shoulder = shoulder)
end

_recall(pk, truth) = count(t -> any(abs.(pk.q .- t) .<= 2HYB_DQ), truth)

@testset "hybrid: default method is unchanged" begin
    golden = Dict(
        "example_tot.dat"     => [110, 138, 160, 199, 234, 248],
        "cubic_tot.dat"       => [94, 123, 138, 153, 172, 222, 268, 306],
        "form-factor_tot.dat" => Int[],
    )
    for (name, idx) in golden
        A = readdlm(joinpath(@__DIR__, "data", name))
        pk = findpeaks(A[:, 1], A[:, 2], A[:, 3])
        @test pk.indices == idx
        @test pk == findpeaks(A[:, 1], A[:, 2], A[:, 3]; method = :persistence)
        @test !haskey(pk, :z)
    end
    @test_throws ArgumentError findpeaks(HYB_Q, ones(HYB_N), ones(HYB_N); method = :nope)
end

@testset "hybrid: synthetic recall and shape" begin
    q, I, σ, t = _hybrid_trace()
    base = findpeaks(q, I, σ)
    hyb  = findpeaks(q, I, σ; method = :hybrid)
    truth = vcat(t.strong, t.weak, t.shoulder)

    @test _recall(hyb, truth) >= _recall(base, truth)
    @test _recall(hyb, t.strong) == length(t.strong)
    @test _recall(hyb, t.weak) == 1
    @test _recall(hyb, t.shoulder) == 1
    @test _recall(base, vcat(t.weak, t.shoulder)) < 2   # the case the hybrid exists for

    @test issorted(hyb.indices)
    @test all(1 .<= hyb.indices .<= HYB_N)
    @test hyb.q == q[hyb.indices]
    @test length(hyb.prominence) == length(hyb.indices) == length(hyb.sharpness)
    @test length(hyb.sharpness_full) == HYB_N
    @test length(hyb.z) == HYB_N
    @test all(hyb.sharpness .> 0)
    @test all(hyb.prominence .>= 0)
    # no spurious peaks beyond a small allowance
    spurious = count(x -> !any(abs(x - t0) <= 3HYB_DQ for t0 in truth), hyb.q)
    @test spurious <= 2
end

@testset "hybrid: few peaks on a trace without any" begin
    for seed in 1:5
        q, I, σ, _ = _hybrid_trace(; seed, peaks = false)
        hyb = findpeaks(q, I, σ; method = :hybrid)
        @test length(hyb.indices) <= 2
        @test all(hyb.sharpness .> 0)
    end
end

@testset "hybrid: deterministic" begin
    q, I, σ, _ = _hybrid_trace()
    a = findpeaks(q, I, σ; method = :hybrid)
    b = findpeaks(q, I, σ; method = :hybrid)
    @test isequal(a, b)
end

@testset "hybrid: keywords" begin
    q, I, σ, _ = _hybrid_trace()
    loose  = findpeaks(q, I, σ; method = :hybrid, k = 8.0)
    strict = findpeaks(q, I, σ; method = :hybrid, k = 1e6, kb = 0.0)
    @test length(loose.indices) >= length(strict.indices)
    # k huge + kb = 0: only persistence peaks remain, none added
    pers = findpeaks(q, I, σ; q_trim_high = 0.07)
    @test strict.indices == pers.indices
    @test findpeaks(q, I; method = :hybrid).indices isa Vector{Int}
end

# A broad peak (half-width ~15 points) must not come back as two peaks.
@testset "hybrid: a broad peak is not reported twice" begin
    q = HYB_Q
    dup = 0
    for seed in 1:200
        rng = MersenneTwister(seed)
        I = 100.0 .+ 400.0 ./ (1 .+ ((q .- 0.2) ./ 0.006) .^ 2)
        σ = sqrt.(I)
        I .+= σ .* randn(rng, HYB_N)
        hyb = findpeaks(q, I, σ; method = :hybrid)
        dup += count(abs.(hyb.q .- 0.2) .< 0.04) > 1
    end
    @test dup == 0
end

@testset "hybrid: added peaks respect q_trim_high" begin
    q, I, σ, _ = _hybrid_trace()
    for trim in (0.07, 0.3)
        hyb = findpeaks(q, I, σ; method = :hybrid, q_trim_high = trim)
        qcut = q[1] + (1 - trim) * (q[end] - q[1])
        @test all(hyb.q .<= qcut)
    end
    # non-uniform grid: the cut is in q, not in points
    qn = q .^ 2 ./ q[end]
    hyb = findpeaks(qn, I, σ; method = :hybrid)
    @test all(hyb.q .<= qn[1] + 0.93 * (qn[end] - qn[1]))
end

@testset "hybrid: persistence peaks near the trace ends are kept" begin
    n = 700
    q = collect(range(0.01, 0.4; length = n))
    I = fill(100.0, n)
    for c in (6, 15, 200, 400, 600)
        I .+= 300.0 ./ (1 .+ ((1:n) .- c) .^ 2 ./ 4)
    end
    σ = sqrt.(I)
    I .+= σ .* randn(MersenneTwister(1), n)
    base = findpeaks(q, I, σ)
    hyb = findpeaks(q, I, σ; method = :hybrid)
    @test 6 in base.indices
    @test all(i -> i in hyb.indices, base.indices)

    # short traces: the widest kernels do not fit; must not error or lose the peak
    qs = q[1:120]; Is = I[1:120]; σs = σ[1:120]
    short = findpeaks(qs, Is, σs; method = :hybrid)
    @test all(i -> i in short.indices, findpeaks(qs, Is, σs).indices)
    @test all(isfinite, short.sharpness)
end

@testset "hybrid: nonpositive intensities do not create peaks" begin
    q, I0, σ0, _ = _hybrid_trace(; seed = 3, peaks = false)
    @test isempty(findpeaks(q, I0, σ0; method = :hybrid).indices)

    beam = copy(I0); beam[1:15] .= 0.0
    @test isempty(findpeaks(q, beam, σ0; method = :hybrid).indices)

    iso = copy(I0); iso[[100, 400, 401, 700]] .= 0.0
    @test isempty(findpeaks(q, iso, σ0; method = :hybrid).indices)

    neg = I0 .- (I0[end] + 4 * σ0[end])      # high-q tail goes slightly negative
    @test isempty(findpeaks(q, neg, σ0; method = :hybrid).indices)

    nanpt = copy(I0); nanpt[50] = NaN
    @test isempty(findpeaks(q, nanpt, σ0; method = :hybrid).indices)
end

@testset "hybrid: result does not depend on intensity scale" begin
    q, I, σ, _ = _hybrid_trace()
    a = findpeaks(q, I, σ; method = :hybrid)
    b = findpeaks(q, I .* 1e-9, σ .* 1e-9; method = :hybrid)
    @test a.indices == b.indices
    c = findpeaks(q, I .* 1e6, σ .* 1e6; method = :hybrid)
    @test a.indices == c.indices
end

@testset "hybrid: without σ, an estimate is used" begin
    q, I, σ, t = _hybrid_trace()
    with_σ = findpeaks(q, I, σ; method = :hybrid)
    no_σ = findpeaks(q, I; method = :hybrid)
    @test count(<=(0.03), no_σ.q) == 0
    @test length(no_σ.indices) <= length(with_σ.indices) + 1
    @test _recall(no_σ, t.strong) == length(t.strong)
end

@testset "hybrid: nonpositive sharpness is replaced by the median of the others" begin
    # A broad shoulder on the convex low-q flank has nonpositive SG sharpness.
    hit = 0
    for seed in 1:8
        q, I, σ, _ = _hybrid_trace(; seed, peaks = false,
            extra = [(0.06, 900.0, 0.008), (0.115, 400.0, 0.0012), (0.16, 400.0, 0.0012),
                     (0.25, 400.0, 0.0012)])
        hyb = findpeaks(q, I, σ; method = :hybrid)
        sf = hyb.sharpness_full[hyb.indices]
        sub = findall(<=(0), sf)
        keep = setdiff(eachindex(sf), sub)
        @test all(hyb.sharpness .> 0)
        @test hyb.sharpness[keep] == sf[keep]
        if !isempty(sub)
            hit += 1
            @test all(hyb.sharpness[sub] .== median(sf[keep]))
        end
    end
    @test hit >= 1                            # the substitution branch is exercised
end

@testset "hybrid: real fixtures" begin
    floors = Dict("example_tot.dat" => 6, "cubic_tot.dat" => 8, "form-factor_tot.dat" => 0)
    spurious_ceiling = Dict("example_tot.dat" => 2, "cubic_tot.dat" => 3, "form-factor_tot.dat" => 2)
    for name in keys(LABELED_PEAKS)
        A = readdlm(joinpath(@__DIR__, "data", name))
        q, I, σ = A[:, 1], A[:, 2], A[:, 3]
        tol = 2 * median(diff(q))
        expected = LABELED_PEAKS[name]
        base = findpeaks(q, I, σ)
        hyb  = findpeaks(q, I, σ; method = :hybrid)
        rec(pk) = count(l -> any(abs.(pk.q .- l) .< tol), expected)
        @test rec(hyb) >= rec(base)
        @test rec(hyb) >= floors[name]
        spur = count(x -> !any(abs.(x .- expected) .< tol), hyb.q)
        @test spur <= spurious_ceiling[name]
        @test issorted(hyb.indices)
        @test all(hyb.sharpness .> 0)
    end
end
