using Test, Random, Statistics, DBInterface, Tables, SQLite
using PeakClassifier, Himalaya
import HimalayaUI

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

function synth_exposure(rng, id, group)
    peaks = [(0.08 + 0.03 * rand(rng) + 0.07 * (k - 1) + 0.02 * k, 0.5 + 2rand(rng), 0.0012) for k in 1:5]
    q, I, σ = synth_trace(rng, peaks)
    LabeledExposure(id, group, q, I, σ, [p[1] for p in peaks], Float64[], Float64[], true, false), peaks
end

@testset "features on a hand-built trace" begin
    rng = MersenneTwister(1)
    q, I, σ = synth_trace(rng, [(0.2, 3.0, 0.002)])
    idx, X = PeakClassifier.candidate_features(q, I, σ)
    @test size(X) == (length(idx), length(FEATURE_NAMES))
    @test all(isfinite, X)
    near = argmin(abs.(q[idx] .- 0.2))
    @test abs(q[idx[near]] - 0.2) < 0.005
    zb = X[:, findfirst(==("z_best"), FEATURE_NAMES)]
    @test zb[near] == maximum(zb)
    @test X[near, findfirst(==("rel_height"), FEATURE_NAMES)] > 0.5
end

@testset "short traces yield no candidates" begin
    q = collect(range(0.1, 0.5; length = 50))
    idx, X = PeakClassifier.candidate_features(q, fill(100.0, 50), ones(50))
    @test isempty(idx) && size(X, 1) == 0
end

@testset "candidate labels" begin
    q = collect(range(0.0, 1.0; length = 1001))
    idx = [100, 300, 500, 700]
    lab(; kw...) = PeakClassifier.label_candidates(q, idx, get(kw, :adds, Float64[]),
        get(kw, :kept, Float64[]), get(kw, :excluded, Float64[]), get(kw, :curated, true), get(kw, :ff, false))
    @test lab(adds = [q[102]], kept = [q[300]], excluded = [q[500]]) == [1, 1, 0, -1]
    @test lab(curated = false) == [-2, -2, -2, -2]
    @test lab(ff = true, curated = false) == [0, 0, 0, 0]
    # a second candidate inside the tolerance of a positive is dropped
    @test lab(adds = [q[100]], kept = Float64[])[1:2] == [1, -1]
    idx2 = [100, 103]
    @test PeakClassifier.label_candidates(q, idx2, [q[100]], Float64[], Float64[], true, false) == [1, -2]
end

@testset "end to end on synthetic traces" begin
    rng = MersenneTwister(7)
    exps = LabeledExposure[]; truth = Dict{Int, Vector{Float64}}()
    for g in 1:4, k in 1:8
        e, pk = synth_exposure(rng, 10g + k, g)
        push!(exps, e); truth[e.id] = [p[1] for p in pk]
    end
    ds = build_dataset(exps)
    @test count(==(1), ds.label) > 0 && count(==(0), ds.label) + count(==(-1), ds.label) > 0

    r = cross_validate(exps, ds; verbose = false, n_trees = 60, thresholds = [0.5])
    @test r.candidate_auc > 0.9
    c = r.by_threshold[1][2]
    @test c.add_recall > 0.8
    @test c.extra_per_curated < r.candidates.extra_per_curated

    held = filter(e -> e.group == 4, exps)
    train_exps = filter(e -> e.group != 4, exps)
    m = train(build_dataset(train_exps); n_trees = 60)
    path = joinpath(mktempdir(), "m.jls")
    save_model(path, m)
    m2 = load_model(path)
    hits = 0; total = 0; extras = 0
    for e in held
        got = predict_peaks(m2, e.q, e.I, e.σ)
        @test got == predict_peaks(m, e.q, e.I, e.σ)
        for t in truth[e.id]
            total += 1; hits += any(abs.(got .- t) .< 0.004)
        end
        extras += count(g -> all(abs.(truth[e.id] .- g) .> 0.004), got)
    end
    @test hits / total > 0.8
    @test extras / length(held) < 2
end

@testset "database loader" begin
    tmp = mktempdir()
    adir = joinpath(tmp, "analysis"); mkpath(adir)
    rng = MersenneTwister(3)
    q, I, σ = synth_trace(rng, [(0.1, 3.0, 0.002), (0.2, 3.0, 0.002), (0.3, 3.0, 0.002)])
    open(joinpath(adir, "ex_tot.dat"), "w") do io
        for i in eachindex(q); println(io, q[i], " ", I[i], " ", σ[i]); end
    end
    db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
    exp_id = HimalayaUI.init_experiment!(db; path = tmp, data_dir = joinpath(tmp, "data"), analysis_dir = adir)
    s_id = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "S")
    e_id = HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = s_id, filename = "ex_tot")
    HimalayaUI.analyze_exposure!(db, e_id, adir)
    auto = [r.q for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT q FROM auto_peaks WHERE exposure_id = ? ORDER BY q", [e_id]))]
    @test length(auto) >= 2
    DBInterface.execute(db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', ?)", [e_id, auto[1]])
    DBInterface.execute(db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', ?)", [e_id, 0.4])
    path = db.file
    DBInterface.close!(db)

    got = load_labeled_exposures(path)
    @test length(got) == 1
    e = got[1]
    @test e.curated && !e.form_factor && e.group == exp_id
    @test e.excluded == [auto[1]]
    @test e.adds == [0.4]
    @test sort(e.kept) == auto[2:end]
end
