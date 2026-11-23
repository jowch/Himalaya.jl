using Test, Random, Statistics
using HimalayaPeakModel, Himalaya
const HPM = HimalayaPeakModel

# Lorentzian peaks on a power-law background with Poisson-like noise.
function synth_trace(rng, peaks; n = 700, amp = 400.0)
    q = collect(range(0.02, 0.5; length = n))
    bg = @. 3000 * (q / 0.02)^-1.8 + 40
    I = copy(bg)
    for (q0, a, w) in peaks
        @. I += a * amp * (w^2 / ((q - q0)^2 + w^2))
    end
    σ = sqrt.(I)
    q, I .+ σ .* randn(rng, n), σ
end

# Cycle through label kinds so every label path is exercised:
# 1 adds only, 2 kept + adds, 3 adds + an excluded stored peak, 0 form factor.
function synth_exposure(rng, id, group)
    kind = id % 4
    kind == 0 && begin
        q, I, σ = synth_trace(rng, [])
        return LabeledExposure(id, group, q, I, σ, Float64[], Float64[], Float64[], false, true), []
    end
    peaks = [(0.08 + 0.03 * rand(rng) + 0.07 * (k - 1) + 0.02 * k, 0.5 + 2rand(rng), 0.0012) for k in 1:5]
    q, I, σ = synth_trace(rng, peaks)
    pq = [p[1] for p in peaks]
    adds, kept, excl = kind == 2 ? (pq[4:5], pq[1:3], Float64[]) :
                       kind == 3 ? (pq, Float64[], [0.33]) : (pq, Float64[], Float64[])
    LabeledExposure(id, group, q, I, σ, adds, kept, excl, true, false), peaks
end

@testset "features on a hand-built trace" begin
    rng = MersenneTwister(1)
    q, I, σ = synth_trace(rng, [(0.2, 3.0, 0.002)])
    idx, X = HPM.candidate_features(q, I, σ)
    @test size(X) == (length(idx), length(FEATURE_NAMES))
    @test all(isfinite, X)
    near = argmin(abs.(q[idx] .- 0.2))
    @test abs(q[idx[near]] - 0.2) < 0.005
    zb = X[:, findfirst(==("z_best"), FEATURE_NAMES)]
    @test zb[near] == maximum(zb)
    @test X[near, findfirst(==("rel_height"), FEATURE_NAMES)] > 0.5
end

@testset "non-finite samples give finite features" begin
    rng = MersenneTwister(2)
    q, I, σ = synth_trace(rng, [(0.2, 3.0, 0.002), (0.35, 3.0, 0.002)])
    j = argmin(abs.(q .- 0.2))
    I[j + 1] = NaN; σ[j - 1] = NaN; σ[j] = 0.0; I[argmin(abs.(q .- 0.35)) - 2] = Inf
    idx, X, _ = HPM.candidate_features(q, I, σ)
    @test !isempty(idx)
    @test all(isfinite, X)
end

@testset "short traces yield no candidates" begin
    q = collect(range(0.1, 0.5; length = 50))
    idx, X = HPM.candidate_features(q, fill(100.0, 50), ones(50))
    @test isempty(idx) && size(X, 1) == 0
end

@testset "candidate labels" begin
    q = collect(range(0.0, 1.0; length = 1001))
    idx = [100, 300, 500, 700]
    lab(; kw...) = HPM.label_candidates(q, idx, get(kw, :adds, Float64[]),
        get(kw, :kept, Float64[]), get(kw, :excluded, Float64[]), get(kw, :curated, true), get(kw, :ff, false))
    @test lab(adds = [q[102]], kept = [q[300]], excluded = [q[500]]) == [1, 1, 0, -1]
    @test lab(curated = false) == [-2, -2, -2, -2]
    @test lab(ff = true, curated = false) == [0, 0, 0, 0]
    # a second candidate inside the tolerance of a positive is dropped
    @test lab(adds = [q[100]], kept = Float64[])[1:2] == [1, -1]
    idx2 = [100, 103]
    @test HPM.label_candidates(q, idx2, [q[100]], Float64[], Float64[], true, false) == [1, -2]
end

@testset "candidate labels: peaks only a model found" begin
    q = collect(range(0.0, 1.0; length = 1001))
    idx = [100, 103, 300, 500, 700]
    model_only = [q[300], q[500]]
    # neither the peak nor its neighbours become negatives or positives
    lab = HPM.label_candidates(q, idx, [q[100]], Float64[], Float64[], true, false; kept_model = model_only)
    @test lab == [1, -2, -2, -2, -1]
    # without them the same candidates are weak negatives
    @test HPM.label_candidates(q, idx, [q[100]], Float64[], Float64[], true, false)[3:4] == [-1, -1]
    # an exclude still beats the weak rule but a model-only peak comes first
    lab2 = HPM.label_candidates(q, idx, Float64[], Float64[], [q[300], q[700]], true, false; kept_model = [q[300]])
    @test lab2 == [-1, -1, -2, -1, 0]
    # a hybrid or both peak is a positive
    @test HPM.label_candidates(q, idx, Float64[], [q[300]], Float64[], true, false; kept_model = [q[500]])[3:4] == [1, -2]
    # form factor: model-only peaks are ignored and everything remaining is a strong negative
    @test HPM.label_candidates(q, idx, Float64[], Float64[], Float64[], false, true; kept_model = model_only) == [0, 0, 0, 0, 0]
    e = LabeledExposure(1, 1, q, q, q, Float64[], Float64[], Float64[], true, false)
    @test e.kept_model == Float64[]
end

@testset "end to end on synthetic traces" begin
    rng = MersenneTwister(7)
    exps = LabeledExposure[]; truth = Dict{Int, Vector{Float64}}()
    for g in 1:4, k in 1:8
        e, pk = synth_exposure(rng, 10g + k, g)
        push!(exps, e); truth[e.id] = [p[1] for p in pk]
    end
    ds = build_dataset(exps)
    @test all(l -> count(==(l), ds.label) > 0, (1, 0, -1))
    ff_rows = [r for r in eachindex(ds.label) if exps[ds.exposure[r]].form_factor]
    @test !isempty(ff_rows) && all(r -> ds.label[r] == 0, ff_rows)

    r = cross_validate(exps, ds; verbose = false, n_trees = 60, thresholds = [0.5])
    @test r.candidate_auc > 0.9
    c = r.by_threshold[1][2]
    @test c.add_recall > 0.8
    @test c.extra_per_labelled < r.candidates.extra_per_labelled

    held = filter(e -> e.group == 4, exps)
    train_exps = filter(e -> e.group != 4, exps)
    m = train(build_dataset(train_exps); n_trees = 60)
    hits = 0; total = 0; extras = 0
    for e in held
        got = predict_peaks(m, e.q, e.I, e.σ)
        for t in truth[e.id]
            total += 1; hits += any(abs.(got .- t) .< 0.004)
        end
        extras += count(g -> all(abs.(truth[e.id] .- g) .> 0.004), got)
    end
    @test total > 0 && hits / total > 0.8
    @test extras / length(held) < 2

    # ablation: no detector columns
    m3 = train(build_dataset(train_exps); n_trees = 30, detector_features = false)
    @test !any(in(HPM.DETECTOR_FEATURES), m3.feature_names)
    @test !isempty(predict_peaks(m3, held[1].q, held[1].I, held[1].σ))
end

# Values taken from the tool before the code moved into this package (same
# seeds), so a change to features, labels or the forest shows up here.
@testset "predictions are pinned on a synthetic set" begin
    rng = MersenneTwister(7)
    exps = LabeledExposure[]
    for g in 1:4, k in 1:8
        push!(exps, synth_exposure(rng, 10g + k, g)[1])
    end
    held = filter(e -> e.group == 4, exps)
    m = train(build_dataset(filter(e -> e.group != 4, exps)); n_trees = 60)
    pins = Dict(41 => [0.1237, 0.1972, 0.307, 0.3949, 0.4636],
                42 => [0.1244, 0.1992, 0.2871, 0.3888, 0.465],
                43 => [0.1093, 0.215, 0.3029, 0.3798, 0.4636],
                44 => Float64[])
    for e in held[1:4]
        got = predict_peaks(m, e.q, e.I, e.σ)
        @test length(got) == length(pins[e.id])
        @test all(isapprox.(got, pins[e.id]; atol = 1e-4))
    end
end

@testset "an exclude labels candidates near it even when no auto peak matches it" begin
    q = collect(range(0.0, 1.0; length = 1001))
    idx = [100, 300, 500, 700]
    # nothing is kept or added near q[300]; the exclude sits 2 points from it
    lab = HPM.label_candidates(q, idx, Float64[], Float64[], [q[302]], true, false)
    @test lab == [-1, 0, -1, -1]
    # farther than the match tolerance from every candidate: no label
    far = HPM.label_candidates(q, idx, Float64[], Float64[], [q[400]], true, false)
    @test far == [-1, -1, -1, -1]
    # a neighbour inside the tolerance of the exclude is negative too
    idx2 = [100, 104]
    @test HPM.label_candidates(q, idx2, Float64[], Float64[], [q[103]], true, false) == [0, 0]

    rng = MersenneTwister(5)
    tq, I, σ = synth_trace(rng, [(0.1, 3.0, 0.002), (0.2, 3.0, 0.002), (0.3, 3.0, 0.002)])
    e = LabeledExposure(1, 1, tq, I, σ, [0.1, 0.2], Float64[], [0.3], true, false)
    ds = build_dataset([e])
    near = argmin(abs.(tq[ds.cand[1]] .- 0.3))
    @test ds.label[near] == 0
end
