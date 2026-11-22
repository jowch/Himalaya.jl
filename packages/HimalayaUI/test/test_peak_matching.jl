using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using Himalaya

# Synthetic grid: step H, peaks near q = 0.1..0.3. Base tolerance at q=0.1 is
# 1e-4 (a tenth of a grid point); the grid-aware tolerance is 6 * H.
const PM_H = 1e-3

function _pm_fixture(; q_step = PM_H, auto_qs = Float64[])
    tmp = mktempdir()
    db  = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
    e_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="x")
    DBInterface.execute(db, "UPDATE exposures SET q_step = ? WHERE id = ?", [q_step, e_id])
    for q in auto_qs
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, ?, 1.0)", [e_id, q])
    end
    (; db, e_id)
end

_curate!(db, e_id, kind, q) = DBInterface.execute(db,
    "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, ?, ?)", [e_id, kind, q])

_grid(n; lo = 0.05, h = PM_H) = collect(range(lo, step = h, length = n))

@testset "peak_match_tol" begin
    @test HimalayaUI.PEAK_MATCH_POINTS == 6
    @test HimalayaUI.peak_match_tol(0.1, nothing) == max(1e-6, 0.1 * 0.001)
    @test HimalayaUI.peak_match_tol(0.1, missing) == max(1e-6, 0.1 * 0.001)
    @test HimalayaUI.peak_match_tol(0.1, 1e-3) ≈ 6e-3
    @test HimalayaUI.peak_match_tol(0.1, 1e-9) == max(1e-6, 0.1 * 0.001)
    @test HimalayaUI.peak_match_tol(0.0, nothing) == 1e-6
end

@testset "migration: exposures.q_step added idempotently" begin
    tmp = mktempdir()
    db  = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    cols() = [r.name for r in Tables.rowtable(DBInterface.execute(db, "PRAGMA table_info(exposures)"))]
    @test "q_step" in cols()
    # Legacy DB: drop the column and the sentinel, then migrate twice.
    DBInterface.execute(db, "ALTER TABLE exposures DROP COLUMN q_step")
    DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?",
        [HimalayaUI.MIGRATION_EXPOSURES_Q_STEP])
    @test !("q_step" in cols())
    HimalayaUI.migrate_exposures_q_step!(db)
    HimalayaUI.migrate_exposures_q_step!(db)
    @test "q_step" in cols()
    HimalayaUI.migrate_schema!(db)
    @test count(==("q_step"), cols()) == 1
end

@testset "analyze_exposure! stores q_step from the trace grid" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    dat = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    cp(dat, joinpath(analysis_dir, "example_tot.dat"))
    db     = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot")
    @test HimalayaUI.exposure_q_step(db, e_id) === nothing
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
    q, _, _ = HimalayaUI.load_dat(dat)
    step = HimalayaUI.exposure_q_step(db, e_id)
    @test step ≈ sum(diff(q)) / (length(q) - 1) rtol = 0.5
    @test step > 0
end

@testset "exclude attaches across a 1-5 grid point move" begin
    for k in 1:5
        fx = _pm_fixture(auto_qs = [0.2 + k * PM_H])
        _curate!(fx.db, fx.e_id, "exclude", 0.2)
        eff = HimalayaUI.effective_peaks(fx.db, fx.e_id, _grid(400), zeros(400))
        @test isempty(eff.q)
    end
    # beyond 6 points it no longer attaches
    fx = _pm_fixture(auto_qs = [0.2 + 7 * PM_H])
    _curate!(fx.db, fx.e_id, "exclude", 0.2)
    @test length(HimalayaUI.effective_peaks(fx.db, fx.e_id, _grid(400), zeros(400)).q) == 1
end

@testset "exclude hides only the nearest of two peaks 3 points apart" begin
    fx = _pm_fixture(auto_qs = [0.2, 0.2 + 3PM_H])
    _curate!(fx.db, fx.e_id, "exclude", 0.2 + 0.4PM_H)
    eff = HimalayaUI.effective_peaks(fx.db, fx.e_id, _grid(400), zeros(400))
    @test eff.q == [0.2 + 3PM_H]
    rows = HimalayaUI.get_peaks_for_exposure(fx.db, fx.e_id)
    @test [Bool(r.excluded) for r in rows if r.source == "auto"] == [true, false]
    @test HimalayaUI.hash_peak_set_from_db(fx.db, fx.e_id) == HimalayaUI.hash_peak_set(eff)
end

@testset "add within tolerance suppresses the auto peak; outside does not" begin
    fx = _pm_fixture(auto_qs = [0.2])
    _curate!(fx.db, fx.e_id, "add", 0.2 + 2PM_H)
    eff = HimalayaUI.effective_peaks(fx.db, fx.e_id, _grid(400), zeros(400))
    @test eff.q == [0.2 + 2PM_H]
    @test eff.peak_kind == [:curation]

    fx = _pm_fixture(auto_qs = [0.2])
    _curate!(fx.db, fx.e_id, "add", 0.2 + 8PM_H)
    eff = HimalayaUI.effective_peaks(fx.db, fx.e_id, _grid(400), zeros(400))
    @test eff.q == [0.2, 0.2 + 8PM_H]
    @test eff.peak_kind == [:auto, :curation]
end

@testset "one add suppresses one auto peak, the nearest" begin
    fx = _pm_fixture(auto_qs = [0.2, 0.2 + 3PM_H])
    _curate!(fx.db, fx.e_id, "add", 0.2 + 2PM_H)
    eff = HimalayaUI.effective_peaks(fx.db, fx.e_id, _grid(400), zeros(400))
    @test eff.q == [0.2, 0.2 + 2PM_H]
end

@testset "NULL q_step keeps the old rule" begin
    fx = _pm_fixture(q_step = nothing, auto_qs = [0.2, 0.3])
    _curate!(fx.db, fx.e_id, "exclude", 0.2 + 0.5e-4)   # inside 0.1% of q
    _curate!(fx.db, fx.e_id, "exclude", 0.3 + 2PM_H)    # outside: must not attach
    _curate!(fx.db, fx.e_id, "add", 0.2 + 0.5e-4)       # add next to auto: old rule keeps both
    eff = HimalayaUI.effective_peaks(fx.db, fx.e_id, _grid(400), zeros(400))
    @test eff.q == [0.2 + 0.5e-4, 0.3]
    @test eff.peak_kind == [:curation, :auto]
    rows = HimalayaUI.get_peaks_for_exposure(fx.db, fx.e_id)
    @test [Bool(r.excluded) for r in rows if r.source == "auto"] == [true, false]
    @test count(r -> r.source == "manual", rows) == 1
end

@testset "NULL q_step: add near auto peak is not suppressed (old behaviour)" begin
    fx = _pm_fixture(q_step = nothing, auto_qs = [0.2])
    _curate!(fx.db, fx.e_id, "add", 0.2 + 0.5e-4)
    eff = HimalayaUI.effective_peaks(fx.db, fx.e_id, _grid(400), zeros(400))
    @test length(eff.q) == 2
end

@testset "UI peak list, effective peaks and hash agree" begin
    fx = _pm_fixture(auto_qs = [0.12, 0.2, 0.2 + 3PM_H, 0.31, 0.4])
    _curate!(fx.db, fx.e_id, "exclude", 0.2 + 0.4PM_H)
    _curate!(fx.db, fx.e_id, "exclude", 0.4 + 5PM_H)
    _curate!(fx.db, fx.e_id, "add", 0.31 + 2PM_H)
    _curate!(fx.db, fx.e_id, "add", 0.5)
    q = _grid(600); I = zeros(600)
    eff = HimalayaUI.effective_peaks(fx.db, fx.e_id, q, I)

    rows = HimalayaUI.get_peaks_for_exposure(fx.db, fx.e_id)
    shown = sort([Float64(r.q) for r in rows if !Bool(r.excluded)])
    @test shown == eff.q

    spec = HimalayaUI._effective_peak_rows(fx.db, fx.e_id)
    @test sort([Float64(r.q) for r in spec]) == eff.q

    with_inproc_routes(fx.db) do call
        r = call("GET", "/api/exposures/$(fx.e_id)/peaks")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test sort([Float64(p.q) for p in list if !p.excluded]) == eff.q
    end

    # Hash agrees with the add-free DB-only path once adds are removed.
    DBInterface.execute(fx.db, "DELETE FROM peak_curations WHERE kind = 'add'")
    eff2 = HimalayaUI.effective_peaks(fx.db, fx.e_id, q, I)
    @test HimalayaUI.hash_peak_set_from_db(fx.db, fx.e_id) == HimalayaUI.hash_peak_set(eff2)
end

@testset "peak toggle route uses nearest matching" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
    step = HimalayaUI.exposure_q_step(db, e_id)
    # Replace the detected peaks by two synthetic ones 3 grid points apart.
    DBInterface.execute(db, "DELETE FROM index_peaks")
    DBInterface.execute(db, "DELETE FROM auto_peaks WHERE exposure_id = ?", [e_id])
    for q in (0.2, 0.2 + 3step)
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, ?, 1.0)", [e_id, q])
    end
    ids = [Int(r.id) for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM auto_peaks ORDER BY q"))]
    hdr = ["Content-Type" => "application/json", "X-Username" => "alice"]
    patch(call, id, ex) = call("PATCH", "/api/peaks/auto/$id"; headers = hdr,
                               body = Vector{UInt8}(JSON3.write(Dict(:excluded => ex))))
    with_inproc_routes(db) do call
        for id in ids
            @test patch(call, id, true).status == 200
        end
        n = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM peak_curations WHERE kind='exclude'"))[1].c
        @test n == 2
        @test patch(call, ids[1], false).status == 200
        list = JSON3.read(String(call("GET", "/api/exposures/$e_id/peaks").body))
        @test [p.excluded for p in list] == [false, true]
    end
end

@testset "diff_update keeps auto ids across a small move" begin
    fx = _pm_fixture(auto_qs = [0.2, 0.3])
    before = Dict(Float64(r.q) => Int(r.id) for r in Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT id, q FROM auto_peaks")))
    I = ones(500)
    pr = (q = [0.2 + 4PM_H, 0.3 - PM_H], indices = [10, 20],
          prominence = [1.0, 1.0], sharpness = [1.0, 1.0])
    HimalayaUI.diff_update_auto_peaks!(fx.db, fx.e_id, pr, I; q_step = PM_H)
    after = Dict(Float64(r.q) => Int(r.id) for r in Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT id, q FROM auto_peaks")))
    @test after[0.2 + 4PM_H] == before[0.2]
    @test after[0.3 - PM_H] == before[0.3]

    # with unknown q_step the same move creates new ids (old rule)
    fx2 = _pm_fixture(auto_qs = [0.2])
    old = Int(Tables.rowtable(DBInterface.execute(fx2.db, "SELECT id FROM auto_peaks"))[1].id)
    HimalayaUI.diff_update_auto_peaks!(fx2.db, fx2.e_id,
        (q = [0.2 + 4PM_H], indices = [10], prominence = [1.0], sharpness = [1.0]), I; q_step = nothing)
    ids = [Int(r.id) for r in Tables.rowtable(DBInterface.execute(fx2.db, "SELECT id FROM auto_peaks"))]
    @test !(old in ids)
end

@testset "diff_update matches one-to-one" begin
    fx = _pm_fixture(auto_qs = [0.2])
    old = Int(Tables.rowtable(DBInterface.execute(fx.db, "SELECT id FROM auto_peaks"))[1].id)
    pr = (q = [0.2 + 2PM_H, 0.2 + 3PM_H], indices = [10, 11],
          prominence = [1.0, 1.0], sharpness = [1.0, 1.0])
    HimalayaUI.diff_update_auto_peaks!(fx.db, fx.e_id, pr, ones(500); q_step = PM_H)
    rows = Tables.rowtable(DBInterface.execute(fx.db, "SELECT id, q FROM auto_peaks ORDER BY q"))
    @test length(rows) == 2
    @test Int(rows[1].id) == old   # nearest new detection keeps the id
end

@testset "replay leaves curation q values as recorded" begin
    fx = _pm_fixture(auto_qs = [0.2])
    req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
    HimalayaUI.apply_event!(fx.db, req; kind = "peak_excluded", entity_type = "exposure",
        entity_id = fx.e_id, payload = Dict(:q => 0.2 - 3PM_H, :auto_peak_id => 1))
    HimalayaUI.apply_event!(fx.db, req; kind = "peak_added", entity_type = "exposure",
        entity_id = fx.e_id, payload = Dict(:q => 0.5))
    HimalayaUI.apply_event!(fx.db, req; kind = "peak_excluded", entity_type = "exposure",
        entity_id = fx.e_id, payload = Dict(:q => 0.7, :auto_peak_id => 9))
    HimalayaUI.apply_event!(fx.db, req; kind = "peak_unexcluded", entity_type = "exposure",
        entity_id = fx.e_id, payload = Dict(:q => 0.7, :auto_peak_id => 9))
    snap() = [(r.kind, Float64(r.q)) for r in Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT kind, q FROM peak_curations WHERE exposure_id = ? ORDER BY q", [fx.e_id]))]
    before = snap()
    @test before == [("exclude", 0.2 - 3PM_H), ("add", 0.5)]
    DBInterface.execute(fx.db, "DELETE FROM peak_curations WHERE exposure_id = ?", [fx.e_id])
    HimalayaUI.rebuild_views_from_log!(fx.db, fx.e_id)
    @test snap() == before
end

@testset "peak_unexcluded removes only the nearest exclude" begin
    fx = _pm_fixture(auto_qs = [0.2, 0.2 + 3PM_H])
    req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
    for q in (0.2, 0.2 + 3PM_H)
        HimalayaUI.apply_event!(fx.db, req; kind = "peak_excluded", entity_type = "exposure",
            entity_id = fx.e_id, payload = Dict(:q => q, :auto_peak_id => 1))
    end
    HimalayaUI.apply_event!(fx.db, req; kind = "peak_unexcluded", entity_type = "exposure",
        entity_id = fx.e_id, payload = Dict(:q => 0.2, :auto_peak_id => 1))
    qs = [Float64(r.q) for r in Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT q FROM peak_curations WHERE kind = 'exclude'"))]
    @test qs == [0.2 + 3PM_H]
end

@testset "peak_unexcluded never deletes a neighbour's exclude" begin
    fx = _pm_fixture(auto_qs = [0.2, 0.2 + 3PM_H])
    req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
    # excludes sit at slightly shifted positions; the nearer one belongs to each peak
    _curate!(fx.db, fx.e_id, "exclude", 0.2 - 0.5PM_H)
    _curate!(fx.db, fx.e_id, "exclude", 0.2 + 3.5PM_H)
    HimalayaUI.apply_event!(fx.db, req; kind = "peak_unexcluded", entity_type = "exposure",
        entity_id = fx.e_id, payload = Dict(:q => 0.2 + 3PM_H, :auto_peak_id => 2))
    qs = [Float64(r.q) for r in Tables.rowtable(DBInterface.execute(fx.db,
        "SELECT q FROM peak_curations WHERE kind = 'exclude'"))]
    @test qs == [0.2 - 0.5PM_H]
    rows = HimalayaUI.get_peaks_for_exposure(fx.db, fx.e_id)
    @test [Bool(r.excluded) for r in rows] == [true, false]
end

@testset "unexclude route on a non-excluded peak is a no-op" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
    step = HimalayaUI.exposure_q_step(db, e_id)
    DBInterface.execute(db, "DELETE FROM index_peaks")
    DBInterface.execute(db, "DELETE FROM auto_peaks WHERE exposure_id = ?", [e_id])
    for q in (0.2, 0.2 + 3step)
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, ?, 1.0)", [e_id, q])
    end
    ids = [Int(r.id) for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM auto_peaks ORDER BY q"))]
    _curate!(db, e_id, "exclude", 0.2)          # only the first peak is excluded
    hdr = ["Content-Type" => "application/json", "X-Username" => "alice"]
    n_events() = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM user_actions WHERE action = 'peak_unexcluded'"))[1].c
    with_inproc_routes(db) do call
        r = call("PATCH", "/api/peaks/auto/$(ids[2])"; headers = hdr,
                 body = Vector{UInt8}(JSON3.write(Dict(:excluded => false))))
        @test r.status == 200
        @test n_events() == 0
        n = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM peak_curations WHERE kind='exclude'"))[1].c
        @test n == 1
        list = JSON3.read(String(call("GET", "/api/exposures/$e_id/peaks").body))
        @test [p.excluded for p in list] == [true, false]
    end
end
