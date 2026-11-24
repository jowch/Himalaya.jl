using Test, Random, Statistics
import JSON
import DecisionTree
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

    r = cross_validate(exps, default_recipe(n_trees = 60); ds, verbose = false, thresholds = [0.5])
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
    @test size(HPM.proba(m3, build_dataset(train_exps)), 1) > 0
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

function synth_set(; groups = 1:4, per = 8, seed = 7)
    rng = MersenneTwister(seed)
    exps = LabeledExposure[]
    for g in groups, k in 1:per
        push!(exps, synth_exposure(rng, 10g + k, g)[1])
    end
    exps
end

json_round(d) = JSON.parse(JSON.json(d))

@testset "forest model survives a JSON round trip exactly" begin
    exps = synth_set()
    ds = build_dataset(exps)
    m = train(ds; n_trees = 40)
    X = ds.X
    p0 = predict_proba(m.fitted, X)
    f2 = from_dict(Fitted, json_round(to_dict(m.fitted)))
    @test f2 isa HPM.ForestFitted
    @test predict_proba(f2, X) == p0
    f3 = from_dict(Fitted, JSON.parse(JSON.json(to_dict(m.fitted)), Dict{String, Any}))
    @test predict_proba(f3, X) == p0
    @test to_dict(m.fitted)["kind"] == "random_forest" && to_dict(m.fitted)["version"] == 2

    # per-tree class fractions and the leaf training counts survive too
    ens0 = m.fitted.forest.ensemble; ens2 = f2.forest.ensemble
    classes = m.fitted.forest.classes
    @test length(ens0.trees) == length(ens2.trees)
    for t in 1:5, r in 1:20:size(X, 1)
        @test DecisionTree.apply_tree_proba(ens0.trees[t], X[r, :], classes) ==
              DecisionTree.apply_tree_proba(ens2.trees[t], X[r, :], classes)
    end
    leaves(n::DecisionTree.Node) = vcat(leaves(n.left), leaves(n.right))
    leaves(l::DecisionTree.Leaf) = [l]
    @test all(t -> sort([sort(l.values) for l in leaves(ens0.trees[t])]) ==
                   sort([sort(l.values) for l in leaves(ens2.trees[t])]) &&
                   [l.majority for l in leaves(ens0.trees[t])] == [l.majority for l in leaves(ens2.trees[t])], 1:5)
end

@testset "loading reads data only" begin
    d = Dict{String, Any}("kind" => "no_such_learner", "version" => 1)
    @test_throws ErrorException from_dict(Fitted, d)
    @test_throws ErrorException from_dict(Learner, d)
    @test_throws ErrorException from_dict(FeatureSet, d)
    bad = to_dict(train(build_dataset(synth_set(groups = 1:2)); n_trees = 3).fitted)
    bad["version"] = 99
    @test_throws ErrorException from_dict(Fitted, bad)
    bad["version"] = 1
    @test_throws ErrorException from_dict(Fitted, bad)
end

@testset "recipe id is stable and sensitive" begin
    r = default_recipe()
    @test recipe_id(r) == recipe_id(from_dict(Recipe, json_round(to_dict(r))))
    @test recipe_id(r) == recipe_id(default_recipe())
    @test length(recipe_id(r)) == 16
    @test recipe_id(r) != recipe_id(default_recipe(n_trees = 100))
    @test recipe_id(r) != recipe_id(default_recipe(detector_features = false))
    @test recipe_id(r) != recipe_id(Recipe(DefaultFeatures(), LogisticLearner()))
    @test to_dict(from_dict(Recipe, json_round(to_dict(r)))) == to_dict(r)
end

@testset "trained model round trip" begin
    exps = synth_set()
    m = train(default_recipe(n_trees = 30, detector_features = false), build_dataset(exps))
    m2 = from_dict(TrainedModel, json_round(to_dict(m)))
    @test m2.feature_names == m.feature_names && m2.trained_on == m.trained_on
    @test m2.threshold == m.threshold && m2.nms == m.nms
    e = exps[1]
    @test predict_peaks(m2, e.q, e.I, e.σ) == predict_peaks(m, e.q, e.I, e.σ)
    @test to_dict(m)["recipe_id"] == recipe_id(m.recipe)
    d = to_dict(m); d["feature_names"] = reverse(d["feature_names"])
    @test_throws ErrorException from_dict(TrainedModel, d)
    @test m2.feature_version == feature_version(m.recipe.features) == m.feature_version
    d = to_dict(m); d["feature_version"] = "default-1;himalaya=0.0.0"
    @test from_dict(TrainedModel, d).feature_version == "default-1;himalaya=0.0.0"
end

@testset "a precomputed hybrid result gives the same features and predictions" begin
    exps = synth_set(groups = 1:2)
    e = exps[1]
    fs = DefaultFeatures()
    idx, X, hyb = candidates_and_features(fs, e.q, e.I, e.σ)
    hres = Himalaya.findpeaks(e.q, e.I, e.σ; method = :hybrid)
    idx2, X2, hyb2 = candidates_and_features(fs, e.q, e.I, e.σ; hybrid = hres)
    @test idx2 == idx && X2 == X && hyb2 === hres && hyb.indices == hres.indices
    m = train(build_dataset(exps); n_trees = 20)
    @test predict_peaks(m, e.q, e.I, e.σ; hybrid = hres) == predict_peaks(m, e.q, e.I, e.σ)
    @test feature_names(DefaultFeatures(false)) == filter(!in(HPM.DETECTOR_FEATURES), FEATURE_NAMES)
    @test occursin(string(pkgversion(Himalaya)), feature_version(fs))
    @test feature_version(fs) != feature_version(DefaultFeatures(false))
    @test to_dict(from_dict(FeatureSet, json_round(to_dict(DefaultFeatures(false))))) == to_dict(DefaultFeatures(false))
end

@testset "a second learner plugs in" begin
    exps = synth_set()
    recipe = Recipe(DefaultFeatures(), LogisticLearner())
    r = cross_validate(exps, recipe; verbose = false, thresholds = [0.5])
    @test r.candidate_auc > 0.8
    held = filter(e -> e.group == 4, exps)
    m = train(recipe, build_dataset(filter(e -> e.group != 4, exps)))
    @test m.fitted isa HPM.LogisticFitted
    ds = build_dataset(held)
    p = predict_proba(m.fitted, ds.X)
    @test all(0 .<= p .<= 1)
    @test HPM.auc(p[ds.label .>= 0], ds.label[ds.label .>= 0] .== 1) > 0.8
    m2 = from_dict(TrainedModel, json_round(to_dict(m)))
    @test predict_proba(m2.fitted, ds.X) == p
    @test all(e -> predict_peaks(m2, e.q, e.I, e.σ) == predict_peaks(m, e.q, e.I, e.σ), held)
    @test any(e -> !isempty(predict_peaks(m, e.q, e.I, e.σ)), held)
end

@testset "a learner registered from outside is used end to end" begin
    struct ConstLearner <: Learner; seed::Int; end
    struct ConstFitted <: Fitted; p::Float64; end
    HPM.fit(l::ConstLearner, X, y; rng) = ConstFitted(mean(y))
    HPM.predict_proba(f::ConstFitted, X) = fill(f.p, size(X, 1))
    HPM.to_dict(l::ConstLearner) = Dict{String, Any}("kind" => "const", "version" => 1, "seed" => l.seed)
    HPM.to_dict(f::ConstFitted) = Dict{String, Any}("kind" => "const", "version" => 1, "p" => f.p)
    register_learner!("const", d -> ConstLearner(d["seed"]), d -> ConstFitted(d["p"]))
    @test_throws ErrorException register_learner!("const", d -> nothing, d -> nothing)
    @test_throws ErrorException register_learner!("random_forest", d -> nothing, d -> nothing)
    @test_throws ErrorException register_features!("default", d -> nothing)
    register_learner!("const", d -> ConstLearner(d["seed"]), d -> ConstFitted(d["p"]); replace = true)
    m = train(Recipe(DefaultFeatures(), ConstLearner(1); threshold = 0.0), build_dataset(synth_set(groups = 1:2)))
    m2 = from_dict(TrainedModel, to_dict(m))
    @test m2.fitted isa ConstFitted && m2.fitted.p == m.fitted.p
end

@testset "logistic fit fails loudly on divergence" begin
    ds = build_dataset(synth_set(groups = 1:2))
    @test_throws ErrorException train(Recipe(DefaultFeatures(), LogisticLearner(lr = 1e300, steps = 50)), ds)
end

@testset "feature cache" begin
    exps = synth_set(groups = 1:2)
    mktempdir() do dir
        key(e) = "trace$(e.id)"
        cache = DiskFeatureCache(dir, feature_version(DefaultFeatures()), key)
        ds = build_dataset(exps; cache)
        files = readdir(dir)
        @test length(files) == length(exps)
        ds2 = build_dataset(exps; cache)
        @test ds2.X == ds.X && ds2.cand == ds.cand && ds2.hybrid == ds.hybrid && ds2.label == ds.label
        @test ds.X == Float64.(Float32.(ds.X))

        path = HPM.cache_path(cache, exps[1])
        idx, X, hyb = HPM.read_features(path, length(feature_names(DefaultFeatures())))
        @test idx == ds.cand[1] && hyb == ds.hybrid[1] && X == ds.X[ds.exposure .== 1, :]
        HPM.write_features(path, idx, X .+ 1, hyb)
        @test build_dataset(exps; cache).X[ds.exposure .== 1, :] == Float64.(Float32.(X .+ 1))

        write(path, read(path)[1:end - 3])
        @test HPM.read_features(path, size(X, 2)) === nothing
        @test build_dataset(exps; cache).X == ds.X
        @test HPM.read_features(path, size(X, 2) + 1) === nothing
        write(path, "junk")
        @test HPM.read_features(path, size(X, 2)) === nothing

        other = DiskFeatureCache(dir, "another-version", key)
        build_dataset(exps; cache = other)
        @test length(readdir(dir)) == 2length(exps)
        @test HPM.cache_path(other, exps[1]) != path

        none = DiskFeatureCache(dir, cache.version, e -> nothing)
        @test HPM.cache_path(none, exps[1]) === nothing
        @test build_dataset(exps; cache = none).X == ds.X
        @test length(readdir(dir)) == 2length(exps)

        blocked = joinpath(dir, "not-a-dir"); write(blocked, "x")
        @test build_dataset(exps; cache = DiskFeatureCache(joinpath(blocked, "sub"), "v", key)).X == ds.X

        pruned = HPM.prune!(cache, Set(key(e) for e in exps[1:2]))
        @test pruned == length(exps) - 2 + length(exps)
        @test sort(readdir(dir)) == sort(vcat(basename.(HPM.cache_path.(Ref(cache), exps[1:2])), "not-a-dir"))
    end
end

@testset "out-of-fold probabilities" begin
    exps = synth_set(groups = 1:3)
    recipe = default_recipe(n_trees = 20)
    ds = build_dataset(exps; features = recipe.features)
    oof = oof_proba(recipe, ds)
    @test length(oof) == length(ds.label) && !any(isnan, oof)
    r = cross_validate(exps, recipe; ds, oof, verbose = false)
    r2 = cross_validate(exps, recipe; ds, verbose = false)
    @test r.by_threshold == r2.by_threshold && r.oof == oof
    @test_throws ErrorException oof_proba(recipe, build_dataset(synth_set(groups = 1:1)))
end

@testset "groups without labelled candidates are predicted by one model trained on everything" begin
    exps = synth_set(groups = 1:3)
    e = exps[1]
    extra = [LabeledExposure(900 + g, 90 + g, e.q, e.I, e.σ, Float64[], Float64[], Float64[], false, false) for g in 1:3]
    recipe = default_recipe(n_trees = 20)
    ds = build_dataset([exps; extra]; features = recipe.features)
    oof = oof_proba(recipe, ds)
    rows = findall(>(90), ds.group)
    @test !isempty(rows) && !any(isnan, oof)
    @test oof[rows] == HPM.proba(train(recipe, ds), ds, rows)
    @test oof[findall(<=(3), ds.group)] == oof_proba(recipe, build_dataset(exps; features = recipe.features))
end

@testset "prune tolerates an unwritable directory and removes stale temporary files" begin
    exps = synth_set(groups = 1:2, per = 2)
    key(e) = "trace$(e.id)"
    mktempdir() do dir
        cache = DiskFeatureCache(dir, "v", key)
        build_dataset(exps; cache)
        tmp = joinpath(dir, "trace1-v.f32.123.1.tmp"); write(tmp, "partial")
        fresh = joinpath(dir, "trace2-v.f32.123.1.tmp"); write(fresh, "in progress")
        run(`touch -d "2 hours ago" $tmp`)
        write(joinpath(dir, "notes.txt"), "keep")
        @test HPM.prune!(cache, Set(key(e) for e in exps)) == 1
        @test "notes.txt" in readdir(dir) && basename(tmp) ∉ readdir(dir) && basename(fresh) in readdir(dir)
        rm(fresh)

        stale = DiskFeatureCache(dir, "older", key)
        build_dataset(exps; cache = stale)
        chmod(dir, 0o500)
        try
            @test HPM.prune!(cache, Set(key(e) for e in exps)) == 0
        finally
            chmod(dir, 0o700)
        end
        @test HPM.prune!(cache, Set(key(e) for e in exps)) == length(exps)
    end
end

@testset "features are rounded to Float32 with or without a cache" begin
    exps = synth_set(groups = 1:2, per = 2)
    plain = build_dataset(exps)
    rounded = build_dataset(exps; float32 = true)
    @test rounded.X == Float64.(Float32.(plain.X)) && rounded.X != plain.X
    mktempdir() do dir
        cached = build_dataset(exps; cache = DiskFeatureCache(dir, "v", e -> "t$(e.id)"))
        @test cached.X == rounded.X
        @test build_dataset(exps; cache = DiskFeatureCache(dir, "v", e -> "t$(e.id)"), float32 = false).X == rounded.X
    end
end
