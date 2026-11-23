using Test, SQLite, DBInterface, Tables
using HimalayaUI: init_experiment!, analyze_exposure!, open_db, create_sample!, create_exposure!

function _beam_off_fixture()
    tmp          = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    open(joinpath(analysis_dir, "dark_tot.dat"), "w") do io
        for i in 1:200
            println(io, 0.005 + 0.001i, " ", i % 7 == 0 ? 2.0 : 0.0, " ", 1.0)
        end
    end
    db     = open_db(joinpath(tmp, "himalaya.db"))
    exp_id = init_experiment!(db; path=tmp, data_dir=joinpath(tmp, "data"), analysis_dir=analysis_dir)
    s_id   = create_sample!(db; experiment_id=exp_id, name="D1")
    lit    = create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot")
    dark   = create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="dark_tot")
    (; db, analysis_dir, exp_id, s_id, lit, dark)
end

stored_median(db, id) = only(Tables.rowtable(DBInterface.execute(db,
    "SELECT trace_median FROM exposures WHERE id = ?", [id]))).trace_median

@testset "beam off" begin
    @testset "trace_median is the middle intensity, nothing for an empty trace" begin
        @test HimalayaUI.trace_median([3.0, 1.0, 2.0]) == 2.0
        @test HimalayaUI.trace_median([4.0, 1.0, 3.0, 2.0]) == 2.5
        @test HimalayaUI.trace_median(Float64[]) === nothing
    end

    @testset "is_beam_off reads the detector first and I₀ only without a trace" begin
        nt = HimalayaUI.NO_TRACE
        @test HimalayaUI.is_beam_off(0.0, missing, 1.2)
        @test !HimalayaUI.is_beam_off(12.0, missing, -0.1)
        @test HimalayaUI.is_beam_off(missing, nt, -0.1)
        @test HimalayaUI.is_beam_off(missing, nt, 0.0)
        @test !HimalayaUI.is_beam_off(missing, nt, 1.2)
        @test !HimalayaUI.is_beam_off(missing, nt, missing)
        @test !HimalayaUI.is_beam_off(missing, missing, -0.1)
        @test !HimalayaUI.is_beam_off(nothing, nothing, nothing)
    end

    @testset "analyze_exposure! stores the trace median" begin
        f = _beam_off_fixture()
        analyze_exposure!(f.db, f.lit, f.analysis_dir)
        analyze_exposure!(f.db, f.dark, f.analysis_dir)
        @test stored_median(f.db, f.lit) > 1
        @test stored_median(f.db, f.dark) == 0.0
    end

    @testset "backfill_trace_medians! fills analyzed images once and skips a missing trace" begin
        f = _beam_off_fixture()
        analyze_exposure!(f.db, f.lit, f.analysis_dir)
        analyze_exposure!(f.db, f.dark, f.analysis_dir)
        DBInterface.execute(f.db, "UPDATE exposures SET trace_median = NULL")
        rm(joinpath(f.analysis_dir, "example_tot.dat"))
        @test HimalayaUI.backfill_trace_medians!(f.db) == 1
        @test stored_median(f.db, f.dark) == 0.0
        @test stored_median(f.db, f.lit) === missing
        @test HimalayaUI.backfill_trace_medians!(f.db) == 0
    end

    @testset "deselect_beam_off! leaves beam-off images out of the selection" begin
        f = _beam_off_fixture()
        analyze_exposure!(f.db, f.lit, f.analysis_dir)
        analyze_exposure!(f.db, f.dark, f.analysis_dir)
        @test HimalayaUI.deselect_beam_off!(f.db, [f.lit, f.dark]) == 1
        @test HimalayaUI.selected_image_ids(f.db, f.s_id) == [f.lit]
        @test HimalayaUI.deselect_beam_off!(f.db, Int[]) == 0
    end

    @testset "the first boot after upgrading leaves beam-off images out of untouched selections, once" begin
        f = _beam_off_fixture()
        analyze_exposure!(f.db, f.lit, f.analysis_dir)
        analyze_exposure!(f.db, f.dark, f.analysis_dir)
        custom = create_sample!(f.db; experiment_id=f.exp_id, name="D2")
        mine = create_exposure!(f.db; experiment_id=f.exp_id, sample_id=custom, filename="dark2")
        theirs = create_exposure!(f.db; experiment_id=f.exp_id, sample_id=custom, filename="lit2")
        DBInterface.execute(f.db, "UPDATE exposures SET trace_median = 0.0 WHERE id = ?", [mine])
        DBInterface.execute(f.db, "UPDATE exposures SET deselected = 1 WHERE id = ?", [theirs])
        DBInterface.execute(f.db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_BEAM_OFF_DESELECTED])
        @test HimalayaUI.deselect_beam_off_once!(f.db) == 1
        @test HimalayaUI.selected_image_ids(f.db, f.s_id) == [f.lit]
        @test HimalayaUI.selected_image_ids(f.db, custom) == [mine]
        DBInterface.execute(f.db, "UPDATE exposures SET deselected = 0 WHERE id = ?", [f.dark])
        @test HimalayaUI.deselect_beam_off_once!(f.db) == 0
        @test HimalayaUI.selected_image_ids(f.db, f.s_id) == [f.lit, f.dark]
    end

    @testset "migrate_exposures_trace_median! adds the column to an old exposures table" begin
        f = _beam_off_fixture()
        DBInterface.execute(f.db, "ALTER TABLE exposures DROP COLUMN trace_median")
        DBInterface.execute(f.db, "DELETE FROM schema_migrations WHERE name = ?",
            [HimalayaUI.MIGRATION_EXPOSURES_TRACE_MEDIAN])
        HimalayaUI.migrate_exposures_trace_median!(f.db)
        HimalayaUI.migrate_exposures_trace_median!(f.db)
        @test "trace_median" in HimalayaUI.cols_of(f.db, "exposures")
    end
end
