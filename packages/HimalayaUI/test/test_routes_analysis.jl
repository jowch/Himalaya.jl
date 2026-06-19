using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using Himalaya

@testset "indices routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    with_inproc_routes(db) do call
        # Indices
        r = call("GET", "/api/exposures/$e_id/indices")
        @test r.status == 200
        indices = JSON3.read(String(r.body))
        @test length(indices) >= 1
        @test haskey(indices[1], :peaks)

        # New in Plan 5: each index has predicted_q (basis × normalized phaseratios)
        # and each matched peak has q_observed attached.
        for entry in indices
            @test haskey(entry, :predicted_q)
            @test length(entry.predicted_q) > 0
            @test all(q -> q > 0, entry.predicted_q)
            for p in entry.peaks
                @test haskey(p, :q_observed)
                @test p.q_observed > 0
            end
        end

        # ngc — present on all indices, finite for the three bicontinuous
        # cubic phases (Ia3d, Pn3m, Im3m), null otherwise. Sign matches
        # core Himalaya.ngc — see src/index.jl.
        cubic = ("Ia3d", "Pn3m", "Im3m")
        for entry in indices
            @test haskey(entry, :ngc)
            if String(entry.phase) in cubic
                @test entry.ngc !== nothing
                @test isfinite(entry.ngc)
            else
                @test entry.ngc === nothing
            end
        end

        # D-10: the legacy /groups routes (GET groups, POST/DELETE members) were
        # retired — the active set is now the durable assignment, covered by the
        # /assignment route tests (test_assignments.jl) and the assignment-native
        # idempotency-replay invariant.
    end
end

# ── M2.4 BACKEND: with_idempotency on speculative routes ─────────────────────

@testset "M2.4: POST /speculative is idempotent under retry with same X-Client-Op-Id" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    peaks = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, q FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 2", [e_id]))
    p1 = Int(peaks[1].id)
    p2 = Int(peaks[2].id)

    with_inproc_routes(db) do call
        body = Dict(:phase => "Lamellar",
                    :anchor_peak_id => p1, :anchor_ratio => 1,
                    :additional => [Dict(:ratio_position => 2, :peak_id => p2)])
        op_id = "uuid-m24-spec-1"

        r1 = call("POST", "/api/exposures/$e_id/speculative";
            headers = ["Content-Type"   => "application/json",
                       "X-Username"     => "alice",
                       "X-Client-Op-Id" => op_id],
            body = Vector{UInt8}(JSON3.write(body)))
        @test r1.status == 200
        body1 = String(copy(r1.body))

        n_after_first = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM indices WHERE exposure_id = ? AND kind = 'speculative'",
            [e_id]))).c
        events_after_first = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions WHERE action = 'speculative_created' AND entity_id = ?",
            [e_id]))).c
        @test n_after_first == 1
        @test events_after_first == 1

        # Same op_id → cached body returned, no new index, no new event row.
        r2 = call("POST", "/api/exposures/$e_id/speculative";
            headers = ["Content-Type"   => "application/json",
                       "X-Username"     => "alice",
                       "X-Client-Op-Id" => op_id],
            body = Vector{UInt8}(JSON3.write(body)))
        @test r2.status == 200
        @test String(copy(r2.body)) == body1

        n_after_second = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM indices WHERE exposure_id = ? AND kind = 'speculative'",
            [e_id]))).c
        events_after_second = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions WHERE action = 'speculative_created' AND entity_id = ?",
            [e_id]))).c
        @test n_after_second == n_after_first
        @test events_after_second == events_after_first
    end
end

@testset "M2.4: DELETE /indices/:id is idempotent under retry with same X-Client-Op-Id" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    peaks = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, q FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 2", [e_id]))
    p1 = Int(peaks[1].id)
    p2 = Int(peaks[2].id)

    new_id = HimalayaUI.insert_speculative_index!(db, e_id, Himalaya.Lamellar,
        Dict{Int,Int}(1 => p1, 2 => p2))

    with_inproc_routes(db) do call
        op_id = "uuid-m24-del-1"

        r1 = call("DELETE", "/api/indices/$new_id";
            headers = ["X-Username"     => "alice",
                       "X-Client-Op-Id" => op_id])
        @test r1.status == 200
        body1 = String(copy(r1.body))

        events_after = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions WHERE action = 'speculative_deleted' AND entity_id = ?",
            [e_id]))).c
        @test events_after == 1

        r2 = call("DELETE", "/api/indices/$new_id";
            headers = ["X-Username"     => "alice",
                       "X-Client-Op-Id" => op_id])
        @test r2.status == 200
        @test String(copy(r2.body)) == body1

        events_after2 = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions WHERE action = 'speculative_deleted' AND entity_id = ?",
            [e_id]))).c
        @test events_after2 == events_after
    end
end
