using Test, Random, DBInterface, SQLite, JSON3
using PeakClassifier
import HimalayaUI

function synth_trace(rng, peaks; n = 700, amp = 400.0)
    q = collect(range(0.02, 0.5; length = n))
    I = @. 3000 * (q / 0.02)^-1.8 + 40
    for (q0, a, w) in peaks
        @. I += a * amp * (w^2 / ((q - q0)^2 + w^2))
    end
    σ = sqrt.(I)
    q, I .+ σ .* randn(rng, n), σ
end

function synth_exposure(rng, id, group)
    peaks = [(0.08 + 0.03 * rand(rng) + 0.07 * (k - 1) + 0.02 * k, 0.5 + 2rand(rng), 0.0012) for k in 1:5]
    q, I, σ = synth_trace(rng, peaks)
    LabeledExposure(id, group, q, I, σ, [p[1] for p in peaks], Float64[], [0.33], true, false)
end

@testset "model file round trip" begin
    rng = MersenneTwister(7)
    exps = [synth_exposure(rng, 10g + k, g) for g in 1:2 for k in 1:6]
    m = train(build_dataset(exps); n_trees = 30)
    path = joinpath(mktempdir(), "sub", "m.json")
    save_model(path, m)
    @test JSON3.read(read(path, String))[:file_version] == 3
    m2 = load_model(path)
    e = exps[1]
    @test predict_peaks(m2, e.q, e.I, e.σ) == predict_peaks(m, e.q, e.I, e.σ)

    m3 = train(build_dataset(exps); n_trees = 30, detector_features = false)
    p3 = joinpath(mktempdir(), "m3.json"); save_model(p3, m3)
    @test predict_peaks(load_model(p3), e.q, e.I, e.σ) == predict_peaks(m3, e.q, e.I, e.σ)

    lr = train(Recipe(DefaultFeatures(), LogisticLearner()), build_dataset(exps))
    pl = joinpath(mktempdir(), "lr.json"); save_model(pl, lr)
    @test predict_peaks(load_model(pl), e.q, e.I, e.σ) == predict_peaks(lr, e.q, e.I, e.σ)

    old = joinpath(mktempdir(), "old.json")
    write(old, JSON3.write(Dict("file_version" => 2, "model" => Dict())))
    @test_throws "retrain" load_model(old)

    serialized = joinpath(mktempdir(), "old.jls")
    open(serialized, "w") do io
        write(io, UInt8[0x37, 0x4a, 0x4c, 0x07, 0x00, 0x01, 0x02])
    end
    @test_throws "retrain" load_model(serialized)

    unreadable = joinpath(mktempdir(), "unreadable.json")
    write(unreadable, "not a model")
    @test_throws "retrain" load_model(unreadable)
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
    db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
    exp_id = HimalayaUI.init_experiment!(db; path = tmp, data_dir = joinpath(tmp, "data"), analysis_dir = adir)
    s_id = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "S")
    e_id = HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = s_id, filename = "ex_tot")
    HimalayaUI.analyze_exposure!(db, e_id, adir)
    DBInterface.execute(db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.4)", [e_id])
    HimalayaUI.record_review!(db, e_id; direct = true)
    (; db, path = db.file, e_id)
end

@testset "loader opens the database read-only" begin
    tmp = mktempdir()
    missing_path = joinpath(tmp, "nope.db")
    @test_throws SQLite.SQLiteException load_labeled_exposures(missing_path)
    @test !isfile(missing_path)

    b = build_loader_db(mktempdir())
    DBInterface.close!(b.db)
    got = load_labeled_exposures(b.path; io = devnull)
    @test [e.id for e in got] == [b.e_id] && got[1].adds == [0.4]
    ro = open_readonly(b.path)
    @test_throws SQLite.SQLiteException DBInterface.execute(ro, "DELETE FROM peak_curations")
    DBInterface.close!(ro)

    # WAL database with frames not yet checkpointed: copy main and -wal while
    # a writer is open, then check the loader leaves both untouched.
    src = mktempdir(); dst = mktempdir()
    b = build_loader_db(src)
    DBInterface.close!(b.db)
    w = HimalayaUI.open_db(b.path)
    DBInterface.execute(w, "PRAGMA wal_autocheckpoint = 0")
    DBInterface.execute(w, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.45)", [b.e_id])
    HimalayaUI.record_review!(w, b.e_id; direct = true)
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
