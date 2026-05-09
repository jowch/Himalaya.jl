using Test
using HTTP, JSON3
using SQLite, DBInterface, Tables
using HimalayaUI
using HimalayaUI: open_db, picker_samples
# `picker_samples` is reachable via `using HimalayaUI:` because it's a
# module-level binding in `comparisons.jl` (which is `include`d into the
# HimalayaUI module). No `export` is required — same pattern as
# `recently_used_exposures` referenced from `routes_picker.jl:51`.

# `experiments.path / data_dir / analysis_dir` are NOT NULL in the schema;
# the '/tmp' literals are valid-path sentinels that picker_samples never
# reads — they exist solely to satisfy the constraint.
@testset "picker_samples helper" begin
    @testset "selected exposure resolves" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "h.db"))
            DBInterface.execute(db,
                "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
            DBInterface.execute(db,
                "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'S')")
            DBInterface.execute(db,
                "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'f1', 0)")
            DBInterface.execute(db,
                "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (101, 10, 'f2', 1)")

            rows = picker_samples(db, 1)
            @test length(rows) == 1
            @test rows[1][:indexing_exposure_id] == 101
            @test length(rows[1][:all_exposures]) == 2
            @test [e[:id] for e in rows[1][:all_exposures]] == [100, 101]   # ORDER BY id ASC
            @test rows[1][:all_exposures][2][:selected] === true             # bool, not 1
        end
    end

    @testset "no selected falls back to highest id" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "h.db"))
            DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'S')")
            DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (200, 10, 'a', 0)")
            DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (199, 10, 'b', 0)")
            rows = picker_samples(db, 1)
            @test rows[1][:indexing_exposure_id] == 200
        end
    end

    @testset "single exposure, selected=0 — still resolves" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "h.db"))
            DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'S')")
            DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (300, 10, 'x', 0)")
            rows = picker_samples(db, 1)
            @test rows[1][:indexing_exposure_id] == 300
        end
    end

    @testset "zero-exposure sample — included with null indexing_exposure_id" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "h.db"))
            DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'EmptyS')")
            rows = picker_samples(db, 1)
            @test length(rows) == 1
            @test rows[1][:indexing_exposure_id] === nothing
            @test isempty(rows[1][:all_exposures])
        end
    end

    @testset "unknown experiment id → empty list" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "h.db"))
            @test picker_samples(db, 99) == Dict{Symbol, Any}[]
        end
    end

    @testset "multi-experiment isolation" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "h.db"))
            DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'A', '/tmp/a', '/tmp/a/data', '/tmp/a/analysis')")
            DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (2, 'B', '/tmp/b', '/tmp/b/data', '/tmp/b/analysis')")
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'A1')")
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (20, 2, 'B1')")
            DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'a', 0)")
            DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (999, 20, 'b', 0)")  # bigger id, wrong exp
            a = picker_samples(db, 1)
            @test length(a) == 1
            @test [e[:id] for e in a[1][:all_exposures]] == [100]
            @test a[1][:indexing_exposure_id] == 100   # not 999 — global MAX(id) would have leaked
        end
    end

    @testset "orphan exposure (sample_id NULL) excluded" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "h.db"))
            DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'S')")
            DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'a', 1)")
            DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (200, NULL, 'orphan', 0)")
            rows = picker_samples(db, 1)
            @test length(rows[1][:all_exposures]) == 1
            @test rows[1][:all_exposures][1][:id] == 100
        end
    end

    @testset "NULL name and display_name render as null in JSON" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "h.db"))
            DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id) VALUES (10, 1)")  # name+display_name NULL
            DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'f', 1)")
            rows = picker_samples(db, 1)
            @test rows[1][:sample][:name]         === nothing
            @test rows[1][:sample][:display_name] === nothing
        end
    end

    @testset "defensive multi-selected legacy data" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "h.db"))
            DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'S')")
            DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'a', 1)")
            DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (200, 10, 'b', 1)")
            rows = picker_samples(db, 1)
            @test rows[1][:indexing_exposure_id] == 200   # higher id wins among multiple selected
        end
    end
end

@testset "GET /api/experiments/:eid/picker-samples" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "h.db"))
        DBInterface.execute(db, "INSERT INTO experiments (id, name, path, data_dir, analysis_dir) VALUES (1, 'E', '/tmp', '/tmp/data', '/tmp/analysis')")
        DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (10, 1, 'S')")
        DBInterface.execute(db, "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (100, 10, 'f1', 1)")

        with_test_server(db) do port, base
            r = HTTP.get("$base/api/experiments/1/picker-samples")
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test length(body) == 1
            @test body[1].indexing_exposure_id == 100
            @test body[1].all_exposures[1].selected === true   # JSON-shape: bool, not 1
            @test haskey(body[1], :indexing_exposure_id)        # null vs absent key

            # Sanity: zero-exposure sample produces null (not absent).
            DBInterface.execute(db, "INSERT INTO samples (id, experiment_id, name) VALUES (11, 1, 'Empty')")
            r2 = HTTP.get("$base/api/experiments/1/picker-samples")
            body2 = JSON3.read(String(r2.body))
            empty_row = first(filter(b -> b.sample.id == 11, collect(body2)))
            @test empty_row.indexing_exposure_id === nothing
        end
    end
end
