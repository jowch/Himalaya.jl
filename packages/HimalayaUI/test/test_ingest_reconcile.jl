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

    @testset "a failed analysis is kept until the analysis version changes" begin
        dir, data, ana = _recon_tree()
        _recon_exposure!(data, ana, "A_S0001_0_001"; t = _T0, h = 58.9, dat = false)
        write(joinpath(ana, "A_S0001_0_001.dat"), "not a trace\n")
        db, id = _recon_db(dir, data, ana)
        HimalayaUI.scan_and_group!(db, id)
        bad = first(_rows(db, "SELECT analysis_error, analysis_version FROM exposures WHERE experiment_id = ?", [id]))
        @test !ismissing(bad.analysis_error) && bad.analysis_error != HimalayaUI.NO_TRACE
        @test bad.analysis_version == HimalayaUI.ANALYSIS_VERSION
        @test !HimalayaUI.ingest_has_work(db, id)

        cp(_RECON_DAT, joinpath(ana, "A_S0001_0_001.dat"); force = true)
        DBInterface.execute(db, "UPDATE exposures SET analysis_version = 'old' WHERE experiment_id = ?", [id])
        @test HimalayaUI.ingest_has_work(db, id)
        HimalayaUI.scan_and_group!(db, id)
        ok = first(_rows(db, "SELECT trace_hash, analysis_error FROM exposures WHERE experiment_id = ?", [id]))
        @test !ismissing(ok.trace_hash)
        @test ismissing(ok.analysis_error)
    end

    @testset "boot resumes an ingest that a restart interrupted" begin
        dir, data, ana = _recon_tree()
        _recon_exposure!(data, ana, "A_S0001_0_001"; t = _T0, h = 58.9)
        db, id = _recon_db(dir, data, ana)
        HimalayaUI.scan_and_group!(db, id; analyze = false)
        DBInterface.execute(db, "UPDATE experiments SET ingest_status = 'scanning' WHERE id = ?", [id])

        HimalayaUI.resume_ingests!(db)
        @test first(_rows(db, "SELECT ingest_status FROM experiments WHERE id = ?", [id])).ingest_status == "complete"
        @test _n(db, "SELECT COUNT(*) AS n FROM exposures WHERE experiment_id = ? AND trace_hash IS NOT NULL", [id]) == 1
    end
end
