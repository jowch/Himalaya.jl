using Test
using SQLite, DBInterface, Tables
using Dates
using HimalayaUI

const _RECON_DAT = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")

function _recon_prp(path; t, h)
    open(path, "w") do io
        println(io, "Time this file was written: ", Dates.format(t, "dd u yyyy HH:MM:SS"))
        println(io, "Exposure time=15")
        println(io, "Beam energy=9000 eV")
        println(io, "Horizontal position=$h mm")
        println(io, "Detector=Pilatus 1M")
    end
end

function _recon_tree()
    dir = mktempdir()
    data = joinpath(dir, "data"); ana = joinpath(dir, "analysis")
    mkpath(data); mkpath(ana)
    (dir, data, ana)
end

function _recon_exposure!(data, ana, stem; t, h, prp = true, tif = true, dat = true)
    tif && write(joinpath(data, "$stem.tif"), "fake tif")
    prp && _recon_prp(joinpath(data, "$stem.prp"); t, h)
    dat && cp(_RECON_DAT, joinpath(ana, "$stem.dat"); force = true)
end

function _recon_db(dir, data, ana)
    db = HimalayaUI.open_db(joinpath(mktempdir(), "h.db"))
    id = HimalayaUI.create_experiment!(db; name = "r", path = dir, data_dir = data, analysis_dir = ana)
    (db, id)
end

_rows(db, sql, p = []) = Tables.rowtable(DBInterface.execute(db, sql, p))
_n(db, sql, p = []) = Int(first(_rows(db, sql, p)).n)

const _T0 = DateTime(2026, 4, 26, 23, 0, 0)

@testset "ingest reconcile" begin
    @testset "a _tot.dat that lands after its tif is analyzed by the next pass" begin
        dir, data, ana = _recon_tree()
        _recon_exposure!(data, ana, "A_S0001_0_001"; t = _T0, h = 58.9, dat = false)
        _recon_exposure!(data, ana, "A_S0002_0_001"; t = _T0 + Second(30), h = 58.9)
        db, id = _recon_db(dir, data, ana)

        HimalayaUI.scan_and_group!(db, id)
        late = first(_rows(db, "SELECT id, trace_hash, analysis_error FROM exposures WHERE filename = 'A_S0001_0_001'"))
        @test ismissing(late.trace_hash)
        @test !ismissing(late.analysis_error)
        @test !HimalayaUI.ingest_has_work(db, id)

        cp(_RECON_DAT, joinpath(ana, "A_S0001_0_001.dat"))
        @test HimalayaUI.ingest_has_work(db, id)
        HimalayaUI.scan_and_group!(db, id)
        late = first(_rows(db, "SELECT trace_hash, analysis_error, analysis_version FROM exposures WHERE filename = 'A_S0001_0_001'"))
        @test !ismissing(late.trace_hash)
        @test ismissing(late.analysis_error)
        @test late.analysis_version == HimalayaUI.ANALYSIS_VERSION
        @test !HimalayaUI.ingest_has_work(db, id)
    end

    @testset "a tif without its prp waits instead of being grouped blind" begin
        dir, data, ana = _recon_tree()
        _recon_exposure!(data, ana, "A_S0001_0_001"; t = _T0, h = 58.9)
        _recon_exposure!(data, ana, "B_S0002_0_001"; t = _T0 + Second(30), h = 62.85, prp = false)
        db, id = _recon_db(dir, data, ana)

        HimalayaUI.scan_and_group!(db, id; analyze = false)
        @test _n(db, "SELECT COUNT(*) AS n FROM exposures WHERE experiment_id = ?", [id]) == 1
        @test HimalayaUI.ingest_has_work(db, id)

        _recon_prp(joinpath(data, "B_S0002_0_001.prp"); t = _T0 + Second(30), h = 62.85)
        HimalayaUI.scan_and_group!(db, id; analyze = false)
        rows = _rows(db, "SELECT filename, sample_id, prp_path FROM exposures WHERE experiment_id = ? ORDER BY filename", [id])
        @test length(rows) == 2
        @test all(!ismissing(r.prp_path) for r in rows)
        @test rows[1].sample_id != rows[2].sample_id
    end

    @testset "a scan never creates a sample with no exposures" begin
        dir, data, ana = _recon_tree()
        _recon_exposure!(data, ana, "A_S0001_0_001"; t = _T0, h = 58.9)
        _recon_exposure!(data, ana, "B_S0002_0_001"; t = _T0 + Second(30), h = 62.85, prp = false)
        _recon_exposure!(data, ana, "C_S0003_0_001"; t = _T0 + Second(60), h = 66.8, prp = false)
        db, id = _recon_db(dir, data, ana)
        HimalayaUI.scan_and_group!(db, id; analyze = false)
        @test _n(db, """SELECT COUNT(*) AS n FROM samples s WHERE s.experiment_id = ?
                        AND NOT EXISTS (SELECT 1 FROM exposures e WHERE e.sample_id = s.id)""", [id]) == 0
    end

    @testset "a TIF-only folder still ingests without prps" begin
        dir, data, ana = _recon_tree()
        _recon_exposure!(data, ana, "A_S0001_0_001"; t = _T0, h = 58.9, prp = false)
        db, id = _recon_db(dir, data, ana)
        HimalayaUI.scan_and_group!(db, id; analyze = false)
        @test _n(db, "SELECT COUNT(*) AS n FROM exposures WHERE experiment_id = ?", [id]) == 1
    end

    @testset "an exposure analyzed by an older version is re-analyzed" begin
        dir, data, ana = _recon_tree()
        _recon_exposure!(data, ana, "A_S0001_0_001"; t = _T0, h = 58.9)
        db, id = _recon_db(dir, data, ana)
        HimalayaUI.scan_and_group!(db, id)
        @test !HimalayaUI.ingest_has_work(db, id)
        DBInterface.execute(db, "UPDATE exposures SET analysis_version = 'old' WHERE experiment_id = ?", [id])
        @test HimalayaUI.ingest_has_work(db, id)
        HimalayaUI.scan_and_group!(db, id)
        @test _n(db, "SELECT COUNT(*) AS n FROM exposures WHERE experiment_id = ? AND analysis_version = ?",
                 [id, HimalayaUI.ANALYSIS_VERSION]) == 1
        @test _n(db, "SELECT COUNT(*) AS n FROM user_actions WHERE action = 'analyze_run'") == 2
    end

    @testset "a failed analysis is retried when its trace file changes" begin
        dir, data, ana = _recon_tree()
        _recon_exposure!(data, ana, "A_S0001_0_001"; t = _T0, h = 58.9, dat = false)
        write(joinpath(ana, "A_S0001_0_001.dat"), "not a trace\n")
        db, id = _recon_db(dir, data, ana)
        HimalayaUI.scan_and_group!(db, id)
        bad = first(_rows(db, "SELECT analysis_error, analysis_version, analysis_trace_stat FROM exposures WHERE experiment_id = ?", [id]))
        @test !ismissing(bad.analysis_error) && bad.analysis_error != HimalayaUI.NO_TRACE
        @test bad.analysis_version == HimalayaUI.ANALYSIS_VERSION
        @test !ismissing(bad.analysis_trace_stat)
        @test !HimalayaUI.ingest_has_work(db, id)

        cp(_RECON_DAT, joinpath(ana, "A_S0001_0_001.dat"); force = true)
        @test HimalayaUI.ingest_has_work(db, id)
        HimalayaUI.scan_and_group!(db, id)
        ok = first(_rows(db, "SELECT trace_hash, analysis_error, analysis_trace_stat FROM exposures WHERE experiment_id = ?", [id]))
        @test !ismissing(ok.trace_hash)
        @test ismissing(ok.analysis_error)
        @test ismissing(ok.analysis_trace_stat)
        @test !HimalayaUI.ingest_has_work(db, id)
    end

    @testset "a forced pass retries a failed analysis even when the file looks unchanged" begin
        dir, data, ana = _recon_tree()
        _recon_exposure!(data, ana, "A_S0001_0_001"; t = _T0, h = 58.9, dat = false)
        trace = joinpath(ana, "A_S0001_0_001.dat")
        write(trace, "not a trace\n")
        db, id = _recon_db(dir, data, ana)
        HimalayaUI.scan_and_group!(db, id)
        cp(_RECON_DAT, trace; force = true)
        DBInterface.execute(db, "UPDATE exposures SET analysis_trace_stat = ? WHERE experiment_id = ?",
            [HimalayaUI._stat_signature(trace), id])
        @test !HimalayaUI.ingest_has_work(db, id)

        HimalayaUI.scan_and_group!(db, id)
        @test ismissing(first(_rows(db, "SELECT trace_hash FROM exposures WHERE experiment_id = ?", [id])).trace_hash)

        @test HimalayaUI.run_ingest!(db, id; force = true) === :ok
        HimalayaUI.wait_prewarm(id)
        ok = first(_rows(db, "SELECT trace_hash, analysis_error FROM exposures WHERE experiment_id = ?", [id]))
        @test !ismissing(ok.trace_hash)
        @test ismissing(ok.analysis_error)
    end

    @testset "re-analysis after a version change keeps the phase call, added peaks and sample call" begin
        dir, data, ana = _recon_tree()
        _recon_exposure!(data, ana, "A_S0001_0_001"; t = _T0, h = 58.9)
        db, id = _recon_db(dir, data, ana)
        HimalayaUI.scan_and_group!(db, id)
        e = first(_rows(db, "SELECT id, sample_id FROM exposures WHERE experiment_id = ?", [id]))
        eid, sid = Int(e.id), Int(e.sample_id)

        DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [eid])
        DBInterface.execute(db, """
            INSERT INTO assignment_members (exposure_id, index_id)
            SELECT ?, m.index_id FROM index_group_members m JOIN index_groups g ON g.id = m.group_id
             WHERE g.exposure_id = ? AND g.kind = 'auto' AND g.active = 1""", [eid, eid])
        DBInterface.execute(db, "INSERT INTO peak_curations (exposure_id, kind, q) VALUES (?, 'add', 0.123)", [eid])
        DBInterface.execute(db, "INSERT INTO sample_calls (sample_id, state) VALUES (?, 'indexed')", [sid])
        HimalayaUI.analyze_exposure!(db, eid, ana)

        members_sql = """SELECT i.phase, i.basis FROM assignment_members m JOIN indices i ON i.id = m.index_id
                          WHERE m.exposure_id = ? ORDER BY i.phase, i.basis"""
        before = _rows(db, members_sql, [eid])
        @test !isempty(before)

        DBInterface.execute(db, "UPDATE exposures SET analysis_version = 'old' WHERE id = ?", [eid])
        @test HimalayaUI.ingest_has_work(db, id)
        @test HimalayaUI.run_ingest!(db, id) === :ok
        HimalayaUI.wait_prewarm(id)
        @test first(_rows(db, "SELECT analysis_version FROM exposures WHERE id = ?", [eid])).analysis_version ==
              HimalayaUI.ANALYSIS_VERSION

        after = _rows(db, members_sql, [eid])
        @test length(after) == length(before)
        @test [r.phase for r in after] == [r.phase for r in before]
        @test all(isapprox(a.basis, b.basis; rtol = 1e-9) for (a, b) in zip(after, before))
        @test first(_rows(db, "SELECT state FROM assignments WHERE exposure_id = ?", [eid])).state == "indexed"
        @test _n(db, "SELECT COUNT(*) AS n FROM peak_curations WHERE exposure_id = ? AND kind = 'add' AND q = 0.123", [eid]) == 1
        @test first(_rows(db, "SELECT state FROM sample_calls WHERE sample_id = ?", [sid])).state == "indexed"
    end

    @testset "boot resumes an ingest that a restart interrupted" begin
        dir, data, ana = _recon_tree()
        _recon_exposure!(data, ana, "A_S0001_0_001"; t = _T0, h = 58.9)
        db, id = _recon_db(dir, data, ana)
        HimalayaUI.scan_and_group!(db, id; analyze = false)
        DBInterface.execute(db, "UPDATE experiments SET ingest_status = 'scanning' WHERE id = ?", [id])

        HimalayaUI.resume_ingests!(db)
        HimalayaUI.wait_prewarm(id)
        @test first(_rows(db, "SELECT ingest_status FROM experiments WHERE id = ?", [id])).ingest_status == "complete"
        @test _n(db, "SELECT COUNT(*) AS n FROM exposures WHERE experiment_id = ? AND trace_hash IS NOT NULL", [id]) == 1
    end

    @testset "the problems route lists exposures without a trace and failed analyses" begin
        dir, data, ana = _recon_tree()
        _recon_exposure!(data, ana, "A_S0001_0_001"; t = _T0, h = 58.9, dat = false)
        _recon_exposure!(data, ana, "A_S0002_0_001"; t = _T0 + Second(30), h = 58.9, dat = false)
        write(joinpath(ana, "A_S0002_0_001.dat"), "not a trace\n")
        _recon_exposure!(data, ana, "A_S0003_0_001"; t = _T0 + Second(60), h = 58.9)
        db, id = _recon_db(dir, data, ana)
        HimalayaUI.scan_and_group!(db, id)
        with_inproc_routes(db) do call
            r = call("GET", "/api/experiments/$id/problems")
            @test r.status == 200
            body = HimalayaUI.JSON3.read(r.body)
            @test body.ingest_error === nothing
            @test [(String(p.filename), String(p.problem)) for p in body.exposures] ==
                  [("A_S0001_0_001", "no_trace"), ("A_S0002_0_001", "failed")]
            @test all(p -> p.sample_id !== nothing && !isempty(p.reason), body.exposures)
            @test call("GET", "/api/experiments/999999/problems").status == 404
        end
    end

    @testset "a failed ingest pass keeps its error until a pass succeeds" begin
        dir, data, ana = _recon_tree()
        _recon_exposure!(data, ana, "A_S0001_0_001"; t = _T0, h = 58.9)
        db, id = _recon_db(dir, data, ana)
        DBInterface.execute(db, "UPDATE experiments SET image_pattern = 'broken.tif' WHERE id = ?", [id])
        @test HimalayaUI.run_ingest!(db, id) === :failed
        row = first(_rows(db, "SELECT ingest_status, ingest_error FROM experiments WHERE id = ?", [id]))
        @test row.ingest_status == "failed"
        @test !ismissing(row.ingest_error) && !isempty(row.ingest_error)
        with_inproc_routes(db) do call
            body = HimalayaUI.JSON3.read(call("GET", "/api/experiments/$id/problems").body)
            @test body.ingest_error == row.ingest_error
        end
        DBInterface.execute(db, "UPDATE experiments SET image_pattern = NULL WHERE id = ?", [id])
        @test HimalayaUI.run_ingest!(db, id) === :ok
        HimalayaUI.wait_prewarm(id)
        row = first(_rows(db, "SELECT ingest_status, ingest_error FROM experiments WHERE id = ?", [id]))
        @test row.ingest_status == "complete"
        @test ismissing(row.ingest_error)
    end

    @testset "opening a database from before ingest_error adds the column" begin
        path = joinpath(mktempdir(), "old.db")
        db = HimalayaUI.open_db(path)
        DBInterface.execute(db, "ALTER TABLE experiments DROP COLUMN ingest_error")
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_EXPERIMENTS_INGEST_ERROR])
        DBInterface.close!(db)
        db = HimalayaUI.open_db(path)
        cols = [r.name for r in _rows(db, "PRAGMA table_info(experiments)", [])]
        @test "ingest_error" in cols
        DBInterface.close!(db)
        db = HimalayaUI.open_db(path)
        @test count(==("ingest_error"), [r.name for r in _rows(db, "PRAGMA table_info(experiments)", [])]) == 1
    end
end
