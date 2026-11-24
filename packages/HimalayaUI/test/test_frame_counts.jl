using Test, Dates, SQLite, DBInterface, Tables, JSON3, HimalayaUI
isdefined(@__MODULE__, :_as_write_folder!) || include("test_setup_mode.jl")

_fc_frame_line(i; rejected = false) =
    "frame #$i: 0.0013  -38.1  1800.0   -nan" * (rejected ? "  Frame rejected!!" : "") * "\n"

function _fc_sample_log(; total = 10, rejected = 1, summary = true)
    io = IOBuffer()
    print(io, "Variance of frames relative to the first frame for buffer files:\n")
    foreach(i -> print(io, _fc_frame_line(i)), 1:3)
    print(io, "Average variance for buffer:\n1790.9\n\nVariance used for rejection: 1790.9\nFactor used for rejection: 1.5\n")
    print(io, "Variance of frames relative to the first frame for sample files:\n")
    foreach(i -> print(io, _fc_frame_line(i; rejected = i <= rejected)), 1:total)
    summary && print(io, "Buffer frames:10 Rejected:0\nSample frames:$total Rejected:$rejected\n")
    String(take!(io))
end

function _fc_buffer_log(; total = 4, rejected = 1)
    "Variance of frames relative to the first frame for buffer files:\n" *
        join((_fc_frame_line(i; rejected = i <= rejected) for i in 1:total))
end

_fc_log_path(root, stem, suffix = "_0_001") =
    joinpath(root, "analysis", "automatic_analysis", "sastool", "$(stem)$(suffix).log")

_fc_touch(path, age_seconds) =
    run(`touch -t $(Dates.format(now() - Second(age_seconds), "yyyymmddHHMM.SS")) $path`)

function _fc_write_log(root, stem, text; suffix = "_0_001", age_seconds = 3600)
    path = _fc_log_path(root, stem, suffix)
    mkpath(dirname(path))
    write(path, text)
    _fc_touch(path, age_seconds)
    path
end

function _fc_counts(db, sid)
    Dict(String(r.filename) => (ismissing(r.frames_total) ? nothing : (Int(r.frames_total), Int(r.frames_rejected)))
         for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT filename, frames_total, frames_rejected FROM exposures WHERE setup_id = ?", [sid])))
end

function _fc_scanned(; logs = Dict{String,String}())
    root = _as_write_folder!(mktempdir(), vcat(_AS_RUN_A, _AS_RUN_B); at = [("wat", 100), ("pep", 103)])
    for (stem, text) in logs
        _fc_write_log(root, stem, text)
    end
    db, eid, sid = _as_live_row_db(root)
    HimalayaUI.scan_and_group!(db, eid)
    (; root, db, eid, sid)
end

@testset "frame counts of an autosampler acquisition" begin
    @testset "the parser reads the sample summary line" begin
        r = HimalayaUI.parse_frame_log(_fc_sample_log(total = 10, rejected = 1))
        @test (r.total, r.rejected, r.final) == (10, 1, true)
        r = HimalayaUI.parse_frame_log(_fc_sample_log(total = 12, rejected = 0))
        @test (r.total, r.rejected) == (12, 0)
    end

    @testset "the parser counts the frames of a first-buffer log that has no summary" begin
        r = HimalayaUI.parse_frame_log(_fc_buffer_log(total = 4, rejected = 1))
        @test (r.total, r.rejected, r.final) == (4, 1, false)
    end

    @testset "the parser gives nothing for a log it cannot trust" begin
        @test HimalayaUI.parse_frame_log("") === nothing
        @test HimalayaUI.parse_frame_log("not a log\n") === nothing
        @test HimalayaUI.parse_frame_log(_fc_sample_log(summary = false)) === nothing
        @test HimalayaUI.parse_frame_log("Sample frames:0 Rejected:0\n") === nothing
        @test HimalayaUI.parse_frame_log("Sample frames:3 Rejected:4\n") === nothing
        @test HimalayaUI.parse_frame_log("Variance of frames relative to the first frame for buffer files:\n") === nothing
    end

    @testset "a scan fills the counts from the log and leaves a missing log NULL" begin
        (; db, sid) = _fc_scanned(logs = Dict(
            "wat_01B_S101" => _fc_sample_log(total = 10, rejected = 1),
            "wat_01C_S102" => _fc_sample_log(total = 10, rejected = 0),
            "wat_01A_B100" => _fc_buffer_log(total = 4, rejected = 1)))
        c = _fc_counts(db, sid)
        @test c["wat_01B_S101"] == (10, 1)
        @test c["wat_01C_S102"] == (10, 0)
        @test c["wat_01A_B100"] == (4, 1)
        @test c["pep_01A_B103"] === nothing
        @test c["pep_01B_S104"] === nothing
    end

    @testset "rows ingested before the columns existed are filled from the logs on disk" begin
        (; root, db, eid, sid) = _fc_scanned()
        @test all(isnothing, values(_fc_counts(db, sid)))
        _fc_write_log(root, "pep_01B_S104", _fc_sample_log(total = 8, rejected = 3))
        @test HimalayaUI.backfill_frame_counts!(db, eid) == 1
        @test _fc_counts(db, sid)["pep_01B_S104"] == (8, 3)
        @test HimalayaUI.backfill_frame_counts!(db, eid) == 0
    end

    @testset "a stored value is never overwritten" begin
        (; root, db, eid, sid) = _fc_scanned()
        DBInterface.execute(db, "UPDATE exposures SET frames_total = 5, frames_rejected = 0 WHERE filename = 'wat_01B_S101'")
        _fc_write_log(root, "wat_01B_S101", _fc_sample_log(total = 10, rejected = 1))
        HimalayaUI.backfill_frame_counts!(db, eid)
        @test _fc_counts(db, sid)["wat_01B_S101"] == (5, 0)
    end

    @testset "a capillary setup is never touched" begin
        (; root, db, eid, sid) = _fc_scanned()
        DBInterface.execute(db, "UPDATE setups SET mode = 'capillary' WHERE id = ?", [sid])
        _fc_write_log(root, "wat_01B_S101", _fc_sample_log())
        @test HimalayaUI.backfill_frame_counts!(db, eid) == 0
        @test all(isnothing, values(_fc_counts(db, sid)))
    end

    @testset "a first-buffer log still being written waits until it has settled" begin
        (; root, db, eid, sid) = _fc_scanned()
        path = _fc_write_log(root, "wat_01A_B100", _fc_buffer_log(total = 2, rejected = 0); age_seconds = 2)
        @test HimalayaUI.backfill_frame_counts!(db, eid) == 0
        @test _fc_counts(db, sid)["wat_01A_B100"] === nothing
        _fc_touch(path, 3600)
        @test HimalayaUI.backfill_frame_counts!(db, eid) == 1
        @test _fc_counts(db, sid)["wat_01A_B100"] == (2, 0)
    end

    @testset "a log with a summary is read at once, however new it is" begin
        (; root, db, eid, sid) = _fc_scanned()
        _fc_write_log(root, "wat_01B_S101", _fc_sample_log(); age_seconds = 1)
        HimalayaUI.backfill_frame_counts!(db, eid)
        @test _fc_counts(db, sid)["wat_01B_S101"] == (10, 1)
    end

    @testset "the log name follows the setup's first-frame suffix" begin
        root = mktempdir()
        _as_add_stem!(root, "wat_01B_S101"; width = 2)
        _as_add_at!(root, "wat", 101)
        _fc_write_log(root, "wat_01B_S101", _fc_sample_log(total = 7, rejected = 2); suffix = "_0_01")
        _fc_write_log(root, "wat_01B_S101", _fc_sample_log(total = 99, rejected = 0); suffix = "_0_001")
        db, eid, sid = _as_live_row_db(root)
        HimalayaUI.scan_and_group!(db, eid)
        @test _fc_counts(db, sid)["wat_01B_S101"] == (7, 2)
    end

    @testset "the columns come last, in the same order in a fresh and a migrated database" begin
        fresh = HimalayaUI.open_db(joinpath(mktempdir(), "f.db"))
        path = joinpath(mktempdir(), "m.db")
        db = HimalayaUI.open_db(path)
        DBInterface.execute(db, "ALTER TABLE exposures DROP COLUMN frames_rejected")
        DBInterface.execute(db, "ALTER TABLE exposures DROP COLUMN frames_total")
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_EXPOSURES_FRAME_COUNTS])
        SQLite.close(db)
        migrated = HimalayaUI.open_db(path)
        @test HimalayaUI.cols_of(migrated, "exposures") == HimalayaUI.cols_of(fresh, "exposures")
        @test HimalayaUI.cols_of(fresh, "exposures")[end-1:end] == ["frames_total", "frames_rejected"]
    end

    @testset "the exposures route carries the counts" begin
        (; db) = _fc_scanned(logs = Dict("wat_01B_S101" => _fc_sample_log(total = 10, rejected = 1)))
        sample_of(stem) = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id FROM exposures WHERE filename = ?", [stem]))).sample_id
        with_inproc_routes(db) do call
            body = only(JSON3.read(String(call("GET", "/api/samples/$(sample_of("wat_01B_S101"))/exposures").body)))
            @test (body.frames_total, body.frames_rejected) == (10, 1)
            body = only(JSON3.read(String(call("GET", "/api/samples/$(sample_of("pep_01B_S104"))/exposures").body)))
            @test body.frames_total === nothing && body.frames_rejected === nothing
        end
    end
end
