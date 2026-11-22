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
    idx, X = PeakClassifier.candidate_features(q, I, σ)
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
    idx, X, _ = PeakClassifier.candidate_features(q, I, σ)
    @test !isempty(idx)
    @test all(isfinite, X)
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
    @test all(l -> count(==(l), ds.label) > 0, (1, 0, -1))
    ff_rows = [r for r in eachindex(ds.label) if exps[ds.exposure[r]].form_factor]
    @test !isempty(ff_rows) && all(r -> ds.label[r] == 0, ff_rows)

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
    @test total > 0 && hits / total > 0.8
    @test extras / length(held) < 2

    # ablation: no detector columns, same score shape, saved and reloaded
    m3 = train(build_dataset(train_exps); n_trees = 30, detector_features = false)
    @test !any(in(PeakClassifier.DETECTOR_FEATURES), m3.feature_names)
    p3 = joinpath(mktempdir(), "m3.jls"); save_model(p3, m3)
    e = held[1]
    @test predict_peaks(load_model(p3), e.q, e.I, e.σ) == predict_peaks(m3, e.q, e.I, e.σ)
end

function write_dat(path, q, I, σ)
    open(path, "w") do io
        for i in eachindex(q); println(io, q[i], " ", I[i], " ", σ[i]); end
    end
end

function build_loader_db(tmp)
    adir = joinpath(tmp, "analysis"); mkpath(adir)
    rng = MersenneTwister(3)
    q, I, σ = synth_trace(rng, [(0.1, 3.0, 0.002), (0.2, 3.0, 0.002), (0.3, 3.0, 0.002)])
    write_dat(joinpath(adir, "ex_tot.dat"), q, I, σ)
    write_dat(joinpath(adir, "ff_tot.dat"), synth_trace(rng, [])...)
    db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
    exp_id = HimalayaUI.init_experiment!(db; path = tmp, data_dir = joinpath(tmp, "data"), analysis_dir = adir)
    s_id = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "S")
    ff_sample = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "F")
    e_id = HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = s_id, filename = "ex_tot")
    ff_id = HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = ff_sample, filename = "ff_tot")
    DBInterface.execute(db, "INSERT OR REPLACE INTO sample_calls (sample_id, state) VALUES (?, 'form_factor')", [ff_sample])
    HimalayaUI.analyze_exposure!(db, e_id, adir)
    auto = [r.q for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT q FROM auto_peaks WHERE exposure_id = ? ORDER BY q", [e_id]))]
    DBInterface.execute(db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'exclude', ?)", [e_id, auto[1]])
    DBInterface.execute(db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', ?)", [e_id, 0.4])
    (db = db, path = db.file, exp_id, s_id, ff_sample, e_id, ff_id, auto, adir)
end

@testset "database loader" begin
    tmp = mktempdir()
    b = build_loader_db(tmp)
    @test length(b.auto) >= 2
    # extra curated exposures that cannot load: no trace file, and no sample
    gone = HimalayaUI.create_exposure!(b.db; experiment_id = b.exp_id, sample_id = b.s_id, filename = "gone_tot")
    orphan = HimalayaUI.create_exposure!(b.db; experiment_id = b.exp_id, sample_id = b.s_id, filename = "orph_tot")
    cp(joinpath(b.adir, "ex_tot.dat"), joinpath(b.adir, "orph_tot.dat"))
    DBInterface.execute(b.db, "UPDATE exposures SET sample_id = NULL WHERE id = ?", [orphan])
    for id in (gone, orphan)
        DBInterface.execute(b.db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.4)", [id])
    end
    DBInterface.close!(b.db)

    log = IOBuffer()
    got = load_labeled_exposures(b.path; io = log)
    @test sort([e.id for e in got]) == sort([b.e_id, b.ff_id, orphan])
    msg = String(take!(log))
    @test occursin("skipped 1 exposures", msg) && occursin("trace file not found", msg)
    e = got[findfirst(e -> e.id == b.e_id, got)]
    @test e.curated && !e.form_factor && e.group == b.exp_id
    @test e.excluded == [b.auto[1]]
    @test e.adds == [0.4]
    @test sort(e.kept) == b.auto[2:end]
    ff = got[findfirst(e -> e.id == b.ff_id, got)]
    @test ff.form_factor && !ff.curated && isempty(ff.adds)
    @test !any(e -> e.id == b.ff_id, load_labeled_exposures(b.path; include_form_factor = false, io = devnull))

    log = IOBuffer()
    by_sample = load_labeled_exposures(b.path; group_by = :sample, io = log)
    @test sort([e.id for e in by_sample]) == sort([b.e_id, b.ff_id])
    @test Dict(e.id => e.group for e in by_sample) == Dict(b.e_id => b.s_id, b.ff_id => b.ff_sample)
    msg = String(take!(log))
    @test occursin("skipped 2 exposures", msg) && occursin("no experiment or sample id", msg)
end

@testset "loader propagates unexpected errors" begin
    tmp = mktempdir()
    b = build_loader_db(tmp)
    DBInterface.execute(b.db, "DROP TABLE auto_peaks")
    DBInterface.close!(b.db)
    @test_throws SQLite.SQLiteException load_labeled_exposures(b.path; io = devnull)
end

@testset "loader opens the database read-only" begin
    tmp = mktempdir()
    missing_path = joinpath(tmp, "nope.db")
    @test_throws SQLite.SQLiteException load_labeled_exposures(missing_path)
    @test !isfile(missing_path)

    # WAL database with frames not yet checkpointed: copy main and -wal while
    # a writer is open, then check the loader leaves both untouched.
    src = mktempdir(); dst = mktempdir()
    b = build_loader_db(src)
    DBInterface.close!(b.db)
    w = HimalayaUI.open_db(b.path)
    DBInterface.execute(w, "PRAGMA wal_autocheckpoint = 0")
    DBInterface.execute(w, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.45)", [b.e_id])
    @test filesize(b.path * "-wal") > 0
    copy = joinpath(dst, "copy.db")
    cp(b.path, copy); cp(b.path * "-wal", copy * "-wal")
    DBInterface.close!(w)
    main0 = read(copy); wal0 = read(copy * "-wal")
    @test !isempty(wal0)
    got = load_labeled_exposures(copy; io = devnull)
    @test any(e -> e.id == b.e_id && 0.45 in e.adds, got)
    @test isfile(copy * "-wal")
    @test read(copy) == main0
    @test read(copy * "-wal") == wal0
end
