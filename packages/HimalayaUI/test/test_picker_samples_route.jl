using Test
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
end
