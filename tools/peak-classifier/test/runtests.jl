using Test, Random, DBInterface, SQLite
import JSON
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
    @test JSON.parse(read(path, String))[:file_version] == 3
    m2 = load_model(path)
    e = exps[1]
    @test predict_peaks(m2, e.q, e.I, e.σ) == predict_peaks(m, e.q, e.I, e.σ)

    m3 = train(build_dataset(exps); n_trees = 30, detector_features = false)
    p3 = joinpath(mktempdir(), "m3.json"); save_model(p3, m3)
    @test predict_peaks(load_model(p3), e.q, e.I, e.σ) == predict_peaks(m3, e.q, e.I, e.σ)

    lr = train(Recipe(DefaultFeatures(), LogisticLearner()), build_dataset(exps))
    pl = joinpath(mktempdir(), "lr.json"); save_model(pl, lr)
    @test predict_peaks(load_model(pl), e.q, e.I, e.σ) == predict_peaks(lr, e.q, e.I, e.σ)

    stale = JSON.parse(read(path, String), Dict{String, Any})
    stale["model"]["feature_version"] = "default-1;detector_features=true;himalaya=0.0.0"
    ps = joinpath(mktempdir(), "stale.json")
    write(ps, JSON.json(stale))
    @test_throws "retrain" load_model(ps)

    old = joinpath(mktempdir(), "old.json")
    write(old, JSON.json(Dict("file_version" => 2, "model" => Dict())))
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

include(joinpath(pkgdir(HimalayaUI), "test", "peak_training_fixture.jl"))
import HimalayaPeakModel as HPM

# 36 labelled exposures in 3 experiments, every one with a confirmed Lamellar call;
# the 12 of the first experiment are re-labelled Hexagonal, which no detector can match.
function eval_fixture(dir)
    b = _pt_db(dir)
    first_x = b.eids[1:12]
    for e in first_x
        DBInterface.execute(b.db, "UPDATE indices SET phase = 'Hexagonal' WHERE exposure_id = ?", [e])
    end
    DBInterface.close!(b.db)
    b
end

@testset "evaluate scores a detector on a database" begin
    mktempdir() do dir
        b = eval_fixture(dir)
        before = read(b.path)
        r = evaluate_database(b.path; detectors = [:hybrid], io = devnull)
        @test read(b.path) == before
        @test r.counts.exposures == 36 && r.counts.labelled == 36 && r.counts.group_by == :experiment
        e = only(r.entries)
        @test e.name == "hybrid" && e.detector == :hybrid && e.threshold === nothing
        @test e.summary.phase == (24, 24) && e.summary.adds == (144, 144) && e.summary.exposures == 36
        @test e.peak_count == round(Int, e.summary.peaks * 36)
        @test phase_table(e) == Dict("Lamellar" => (24, 24))
        @test r.counts.phase_calls.all == 36 && r.counts.phase_calls.reachable == 24 && r.counts.phase_targets == 24
        @test r.counts.phase_calls.per_phase["Hexagonal"] == (reachable = 0, all = 12)
        @test isempty(r.paired)

        every = evaluate_database(b.path; detectors = [:hybrid], phase_targets = :all, io = devnull)
        @test r.phase_targets == :reachable && every.phase_targets == :all
        @test only(every.entries).summary.phase == (24, 36) && every.counts.phase_targets == 36
        @test phase_table(only(every.entries)) == Dict("Lamellar" => (24, 24), "Hexagonal" => (0, 12))
        @test every.counts.phase_calls == r.counts.phase_calls
    end
end

@testset "evaluate compares detectors against the first, per phase" begin
    mktempdir() do dir
        b = eval_fixture(dir)
        r = evaluate_database(b.path; detectors = [:hybrid, :persistence], io = devnull)
        @test [e.name for e in r.entries] == ["hybrid", "persistence"]
        a, p = r.entries
        pair = only(r.paired)
        @test pair.name == "persistence" && pair.against == "hybrid"
        for (phase, (w, l, both, neither)) in pair.phases
            na = phase_table(a)[phase]; np = phase_table(p)[phase]
            @test w + both == np[1] && l + both == na[1]
            @test w + l + both + neither == na[2]
        end
        @test sum(x -> x[1], values(pair.phases)) - sum(x -> x[2], values(pair.phases)) ==
              p.summary.phase[1] - a.summary.phase[1]
        @test !haskey(pair.phases, "Hexagonal")
    end
end

@testset "evaluate scores a recipe out of fold and sweeps thresholds" begin
    mktempdir() do dir
        b = eval_fixture(dir)
        recipe = HPM.default_recipe(n_trees = 15)
        r = evaluate_database(b.path; detectors = [:hybrid, :recipe], recipe, thresholds = :sweep, io = devnull)
        @test [e.name for e in r.entries][1] == "hybrid"
        @test [e.threshold for e in r.entries[2:end]] == collect(HimalayaUI.THRESHOLD_GRID)
        @test all(e -> e.detector == :recipe && e.summary.exposures == 36, r.entries[2:end])
        @test r.entries[2].name == "recipe@0.3"
        @test length(r.paired) == length(HimalayaUI.THRESHOLD_GRID)
        @test r.entries[end].summary.peaks <= r.entries[2].summary.peaks

        one = evaluate_database(b.path; detectors = [:recipe], recipe, thresholds = 0.6, io = devnull)
        @test only(one.entries).name == "recipe@0.6"
        @test isequal(only(one.entries).summary, r.entries[findfirst(e -> e.threshold == 0.6, r.entries)].summary)
        @test_throws ArgumentError evaluate_database(b.path; detectors = [:nope], io = devnull)
    end
end

@testset "oracle scores the user's own peaks as a reference" begin
    mktempdir() do dir
        b = eval_fixture(dir)
        r = evaluate_database(b.path; detectors = [:oracle, :hybrid], phase_targets = :all, io = devnull)
        o = r.entries[1]
        @test o.name == "oracle" && o.detector == :oracle && o.threshold === nothing
        @test o.summary.adds == (144, 144) && o.summary.excluded == (0, 12)
        @test phase_table(o) == Dict("Lamellar" => (24, 24), "Hexagonal" => (0, 12))
        @test o.peak_count >= 144
        @test only(r.paired).phases["Lamellar"][2] == 0
        @test occursin("oracle", sprint(print_evaluation, r))
    end
end

@testset "evaluate checks its options before loading" begin
    mktempdir() do dir
        b = eval_fixture(dir)
        @test_throws ArgumentError evaluate_database(b.path; phase_auto = :typo, io = devnull)
        @test_throws ArgumentError evaluate_database(b.path; phase_targets = :typo, io = devnull)
        @test_logs (:warn, r"--threshold|threshold") evaluate_database(b.path; detectors = [:hybrid], thresholds = 0.7, io = devnull)
        r = evaluate_database(b.path; detectors = [:recipe], recipe = HPM.default_recipe(n_trees = 10),
                              thresholds = :sweep, io = devnull)
        @test r.recipe_id === nothing
        one = evaluate_database(b.path; detectors = [:recipe], recipe = HPM.default_recipe(n_trees = 10),
                                thresholds = 0.6, io = devnull)
        @test one.recipe_id isa String
    end
end

@testset "evaluate prints and writes the same numbers" begin
    mktempdir() do dir
        b = eval_fixture(dir)
        r = evaluate_database(b.path; detectors = [:hybrid, :persistence], io = devnull)
        text = sprint(print_evaluation, r)
        @test occursin("24/24", text) && occursin("144/144", text) && occursin("Hexagonal 0/12", text)
        @test occursin("24 of 36", text)
        @test occursin("persistence", text) && occursin("wins", text)
        path = joinpath(dir, "sub", "out.json")
        write_evaluation(path, r)
        j = JSON.parse(read(path, String))
        @test [x.name for x in j.detectors] == ["hybrid", "persistence"]
        h = j.detectors[1]
        @test h.metrics.phase.correct == 24 && h.metrics.phase.total == 24
        @test j.phase_targets == "reachable" && j.counts.phase_calls.all == 36 && j.counts.phase_calls.reachable == 24
        @test j.counts.phase_calls.per_phase.Hexagonal.reachable == 0 && j.counts.phase_calls.per_phase.Hexagonal.all == 12
        @test h.metrics.adds.correct == 144 && h.metrics.adds.total == 144
        @test h.metrics.peaks_per_exposure.total == 36
        @test !haskey(h.per_phase, :Hexagonal)
        @test j.paired[1].against == "hybrid" && j.counts.exposures == 36
        @test !occursin(dir, read(path, String))
    end
end

@testset "the training gate scores with the same function as evaluate" begin
    mktempdir() do dir
        b = eval_fixture(dir)
        ro = HimalayaUI.open_readonly_db(b.path)
        run = JSON.parse(JSON.json(HimalayaUI.run_peak_training(ro; recipe = HPM.default_recipe(n_trees = 15),
                                                                  io = devnull, peaks_slack = 100, ff_slack = 100)))
        DBInterface.close!(ro)
        h = only(evaluate_database(b.path; detectors = [:hybrid], io = devnull).entries).summary
        g = run.metrics.comparison.hybrid
        @test g.phase == [h.phase...] && g.adds == [h.adds...] && g.excluded == [h.excluded...]
        @test g.peaks == h.peaks && g.exposures == h.exposures

        m = run.metrics.comparison.model
        thr = run.metrics.comparison.threshold
        mine = only(evaluate_database(b.path; detectors = [:recipe], recipe = HPM.default_recipe(n_trees = 15),
                                      thresholds = thr, io = devnull).entries).summary
        @test m.phase == [mine.phase...] && m.adds == [mine.adds...] && m.excluded == [mine.excluded...]
        @test m.peaks == mine.peaks && m.exposures == mine.exposures
    end
end
