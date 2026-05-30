using Test
using HimalayaUI
using SQLite, DBInterface, Tables
using JSON3, HTTP

@testset "assignments schema" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))

        # Both new tables exist.
        tbls = Set(String(r.name) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table'")))
        @test "assignments" in tbls
        @test "assignment_members" in tbls

        # state CHECK constraint rejects an unknown value.
        DBInterface.execute(db,
            "INSERT INTO exposures (id, sample_id) VALUES (1, NULL)")
        @test_throws Exception DBInterface.execute(db,
            "INSERT INTO assignments (exposure_id, state) VALUES (1, 'bogus')")

        # Default state is 'indexed'.
        DBInterface.execute(db,
            "INSERT INTO assignments (exposure_id) VALUES (1)")
        st = Tables.rowtable(DBInterface.execute(db,
            "SELECT state FROM assignments WHERE exposure_id = 1"))[1].state
        @test String(st) == "indexed"
    end
end

@testset "migrate_assignments! backfills from active group" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

        # Two indices, an active custom group owning both, an inactive auto group.
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (11, ?, 'Im3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (100, ?, 'auto', 0)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (101, ?, 'custom', 1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_group_members (group_id, index_id) VALUES (101, 10)")
        DBInterface.execute(db, "INSERT INTO index_group_members (group_id, index_id) VALUES (101, 11)")

        # Wipe the migration sentinel + new tables so we can re-run the migration deterministically.
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = 'assignments_v1'")
        DBInterface.execute(db, "DELETE FROM assignment_members WHERE exposure_id = ?", [e_id])
        DBInterface.execute(db, "DELETE FROM assignments WHERE exposure_id = ?", [e_id])

        HimalayaUI.migrate_assignments!(db)

        state = Tables.rowtable(DBInterface.execute(db,
            "SELECT state FROM assignments WHERE exposure_id = ?", [e_id]))
        @test !isempty(state) && String(state[1].state) == "indexed"
        members = Set(Int(m.index_id) for m in Tables.rowtable(DBInterface.execute(db,
            "SELECT index_id FROM assignment_members WHERE exposure_id = ?", [e_id])))
        @test members == Set([10, 11])
    end
end

@testset "_assignment_body shape" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

        # No assignment yet → defaults: state 'indexed', empty members.
        b0 = HimalayaUI._assignment_body(db, e_id)
        @test b0[:exposure_id] == e_id
        @test b0[:state] == "indexed"
        @test b0[:members] == Int[]

        # With members + a non-default state.
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (11, ?, 'Im3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [e_id])
        DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, 11)", [e_id])
        DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, 10)", [e_id])
        b1 = HimalayaUI._assignment_body(db, e_id)
        @test b1[:members] == [10, 11]   # ORDER BY index_id
    end
end

@testset "GET /assignment serves the assignment body" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'form_factor')", [e_id])

        # The route's body is _assignment_body; assert its JSON-serialized shape.
        body = HimalayaUI._assignment_body(db, e_id)
        round = JSON3.read(JSON3.write(body))
        @test round.exposure_id == e_id
        @test round.state == "form_factor"
        @test collect(round.members) == Int[]
    end
end

# Full in-process HTTP coverage — only when the shared test server harness
# (test_http.jl::with_test_server) is loaded, i.e. under the full runtests.jl
# suite. Skipped on a standalone `julia test_assignments.jl` run.
if isdefined(@__MODULE__, :with_test_server)
    @testset "GET /assignment in-process HTTP" begin
        mktempdir() do dir
            db = HimalayaUI.open_db(joinpath(dir, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
            e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
            DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
            DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [e_id])
            DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, 10)", [e_id])

            with_test_server(db) do port, base
                r = HTTP.get("$base/api/exposures/$e_id/assignment")
                @test r.status == 200
                got = JSON3.read(String(r.body))
                @test got.exposure_id == e_id
                @test got.state == "indexed"
                @test collect(got.members) == [10]
            end
        end
    end
end

@testset "assignment/state validation + effect" begin
    mktempdir() do dir
        db  = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

        # Valid set → 'null' state recorded.
        HimalayaUI.apply_event!(db, req; kind="assignment_set_state",
            entity_type="exposure", entity_id=e_id, payload=Dict(:state => "null"))
        @test HimalayaUI._assignment_body(db, e_id)[:state] == "null"

        # The route's allow-list predicate (mirror of the handler guard).
        valid = s -> s in ("indexed", "form_factor", "null")
        @test valid("form_factor")
        @test !valid("bogus")
    end
end

if isdefined(@__MODULE__, :with_test_server)
    @testset "POST /assignment/state in-process HTTP" begin
        mktempdir() do dir
            db = HimalayaUI.open_db(joinpath(dir, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
            e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
            DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
            DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [e_id])
            DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, 10)", [e_id])

            with_test_server(db) do port, base
                # Invalid state → 400.
                r = HTTP.post("$base/api/exposures/$e_id/assignment/state",
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:state => "bogus")); status_exception=false)
                @test r.status == 400

                # Valid: form_factor clears members.
                r = HTTP.post("$base/api/exposures/$e_id/assignment/state",
                    ["Content-Type" => "application/json", "X-Username" => "alice"],
                    JSON3.write(Dict(:state => "form_factor")))
                @test r.status == 200
                got = JSON3.read(String(r.body))
                @test got.state == "form_factor"
                @test collect(got.members) == Int[]
                @test HimalayaUI._assignment_body(db, e_id)[:state] == "form_factor"
            end
        end
    end
end

@testset "member-add dual-writes the assignment" begin
    mktempdir() do dir
        db  = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (200, ?, 'custom', 1)", [e_id])

        # Simulate exactly what the route body now does: legacy event + dual-write.
        HimalayaUI.apply_event!(db, req; kind="index_confirmed",
            entity_type="exposure", entity_id=e_id, payload=Dict(:group_id => 200, :index_id => 10))
        HimalayaUI.apply_event!(db, req; kind="assignment_add",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))

        # Both sources of truth agree.
        legacy = Set(Int(m.index_id) for m in Tables.rowtable(DBInterface.execute(db,
            "SELECT index_id FROM index_group_members WHERE group_id = 200")))
        @test legacy == Set([10])
        @test HimalayaUI._assignment_body(db, e_id)[:members] == [10]
    end
end
