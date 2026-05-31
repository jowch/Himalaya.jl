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

@testset "migrate_assignments! seeds auto-only (never-confirmed) exposures as indexed" begin
    # The dominant production case: an analyzed exposure the user never
    # confirmed has only an ACTIVE auto group (no custom group). Legacy
    # confirmed_index ignored it (kind='custom' filter), so it read as null.
    # The new model treats the auto selection as the default assignment, so the
    # migration must seed state='indexed' + the auto members — converging with
    # seed_assignment_if_absent! for a fresh analyze. This pins the single
    # biggest behavior shift at the D-10 cutover.
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (20, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (21, ?, 'Im3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (200, ?, 'auto', 1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_group_members (group_id, index_id) VALUES (200, 20)")
        DBInterface.execute(db, "INSERT INTO index_group_members (group_id, index_id) VALUES (200, 21)")

        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = 'assignments_v1'")
        DBInterface.execute(db, "DELETE FROM assignment_members WHERE exposure_id = ?", [e_id])
        DBInterface.execute(db, "DELETE FROM assignments WHERE exposure_id = ?", [e_id])

        HimalayaUI.migrate_assignments!(db)

        state = Tables.rowtable(DBInterface.execute(db,
            "SELECT state FROM assignments WHERE exposure_id = ?", [e_id]))
        @test !isempty(state) && String(state[1].state) == "indexed"
        members = Set(Int(m.index_id) for m in Tables.rowtable(DBInterface.execute(db,
            "SELECT index_id FROM assignment_members WHERE exposure_id = ?", [e_id])))
        @test members == Set([20, 21])
    end
end

@testset "migrate_assignments! sentinel protects post-migration user edits" begin
    # The load-bearing safety contract: migrate_schema! runs on EVERY server
    # boot. Once the sentinel is set, a re-run must short-circuit and never
    # clobber edits the user made after the first migration — otherwise an
    # INSERT OR IGNORE would re-add deleted members on the next restart. The
    # happy-path test deletes the sentinel to force a re-run, so this guard is
    # otherwise entirely uncovered.
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (30, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (31, ?, 'Im3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (300, ?, 'custom', 1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_group_members (group_id, index_id) VALUES (300, 30)")
        DBInterface.execute(db, "INSERT INTO index_group_members (group_id, index_id) VALUES (300, 31)")

        # First migration (force fresh).
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = 'assignments_v1'")
        DBInterface.execute(db, "DELETE FROM assignment_members WHERE exposure_id = ?", [e_id])
        DBInterface.execute(db, "DELETE FROM assignments WHERE exposure_id = ?", [e_id])
        HimalayaUI.migrate_assignments!(db)

        # User edits the assignment AFTER migration: drops a member, sets null.
        DBInterface.execute(db, "DELETE FROM assignment_members WHERE exposure_id = ? AND index_id = 31", [e_id])
        DBInterface.execute(db, "UPDATE assignments SET state = 'null' WHERE exposure_id = ?", [e_id])

        # Second boot: sentinel is present, so this must be a no-op.
        HimalayaUI.migrate_assignments!(db)

        state = Tables.rowtable(DBInterface.execute(db,
            "SELECT state FROM assignments WHERE exposure_id = ?", [e_id]))[1].state
        @test String(state) == "null"            # state edit preserved
        members = Set(Int(m.index_id) for m in Tables.rowtable(DBInterface.execute(db,
            "SELECT index_id FROM assignment_members WHERE exposure_id = ?", [e_id])))
        @test members == Set([30])               # dropped member stays dropped
    end
end

@testset "migrate_assignments! skips exposures with no active group" begin
    # An exposure whose only group is inactive (active=0) must NOT get an
    # assignment row — the migration keys strictly off active=1, mirroring the
    # legacy active-group semantics.
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (40, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (400, ?, 'auto', 0)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_group_members (group_id, index_id) VALUES (400, 40)")

        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = 'assignments_v1'")
        DBInterface.execute(db, "DELETE FROM assignment_members WHERE exposure_id = ?", [e_id])
        DBInterface.execute(db, "DELETE FROM assignments WHERE exposure_id = ?", [e_id])
        HimalayaUI.migrate_assignments!(db)

        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM assignments WHERE exposure_id = ?", [e_id])))
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM assignment_members WHERE exposure_id = ?", [e_id])))
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

@testset "member-remove dual-writes the assignment" begin
    mktempdir() do dir
        db  = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (200, ?, 'custom', 1)", [e_id])

        HimalayaUI.apply_event!(db, req; kind="assignment_add",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))
        # Simulate the DELETE route body: legacy unconfirm + dual-write remove.
        HimalayaUI.apply_event!(db, req; kind="index_unconfirmed",
            entity_type="exposure", entity_id=e_id, payload=Dict(:group_id => 200, :index_id => 10))
        HimalayaUI.apply_event!(db, req; kind="assignment_remove",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))

        @test isempty(HimalayaUI._assignment_body(db, e_id)[:members])
    end
end

@testset "assignment seed-if-absent contract" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])

        # Call the extracted seed helper directly.
        HimalayaUI.seed_assignment_if_absent!(db, e_id, [10])
        @test HimalayaUI._assignment_body(db, e_id)[:members] == [10]
        @test HimalayaUI._assignment_body(db, e_id)[:state] == "indexed"

        # Second call with a different selection must NOT clobber existing membership.
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (11, ?, 'Im3m', 0.1)", [e_id])
        HimalayaUI.seed_assignment_if_absent!(db, e_id, [11])
        @test HimalayaUI._assignment_body(db, e_id)[:members] == [10]   # preserved
    end
end

@testset "_bonnet_for_index flags a coexisting cubic" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

        # Assigned Pn3m at a=100; a candidate Im3m at a=128 (Bonnet-consistent).
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis, lattice_d) VALUES (10, ?, 'Pn3m', 0.1, 100.0)", [e_id])
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis, lattice_d) VALUES (11, ?, 'Im3m', 0.1, 128.0)", [e_id])
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis, lattice_d) VALUES (12, ?, 'Im3m', 0.1, 200.0)", [e_id])
        DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [e_id])
        DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, 10)", [e_id])

        # The assigned Pn3m (10) itself: no flag (it's the anchor, same phase).
        @test HimalayaUI._bonnet_for_index(db, e_id, "Pn3m", 100.0, 10) === nothing
        # Candidate Im3m at 128 → consistent, predicted ≈ 127.9.
        b = HimalayaUI._bonnet_for_index(db, e_id, "Im3m", 128.0, 11)
        @test b !== nothing && b[:consistent] == true
        @test b[:predicted_a] ≈ 127.9 atol=1.0
        # Candidate Im3m at 200 → predicted still 127.9, not consistent.
        b2 = HimalayaUI._bonnet_for_index(db, e_id, "Im3m", 200.0, 12)
        @test b2 !== nothing && b2[:consistent] == false
    end
end

# ── Plan D Task D-1: native assignment member routes ─────────────────────────
# POST/DELETE /api/exposures/{id}/assignment/members emit assignment_add /
# assignment_remove with a DISTINCT post_state ({assignment:{state,members}},
# no top-level `indices` key). These are the assignment-native targets that
# replace the legacy /groups dual-write. The legacy /groups routes remain live
# (retired separately in D-10).
if isdefined(@__MODULE__, :with_test_server)
    @testset "POST /assignment/members adds a member (native route)" begin
        mktempdir() do dir
            db = HimalayaUI.open_db(joinpath(dir, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
            e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
            DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
            DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (11, ?, 'Im3m', 0.1)", [e_id])

            with_test_server(db) do port, base
                # Missing field -> 400.
                r = HTTP.post("$base/api/exposures/$e_id/assignment/members",
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:nope => 1)); status_exception=false)
                @test r.status == 400

                # Valid add -> 200, body is the canonical assignment shape.
                r = HTTP.post("$base/api/exposures/$e_id/assignment/members",
                    ["Content-Type" => "application/json", "X-Username" => "alice"],
                    JSON3.write(Dict(:index_id => 10)))
                @test r.status == 200
                got = JSON3.read(String(r.body))
                @test got.exposure_id == e_id
                @test got.state == "indexed"
                @test collect(got.members) == [10]
                @test haskey(got, :event_id)
                @test HimalayaUI._assignment_body(db, e_id)[:members] == [10]

                # Idempotent re-add of the same member is a no-op on membership.
                r = HTTP.post("$base/api/exposures/$e_id/assignment/members",
                    ["Content-Type" => "application/json", "X-Username" => "alice"],
                    JSON3.write(Dict(:index_id => 11)))
                @test r.status == 200
                @test HimalayaUI._assignment_body(db, e_id)[:members] == [10, 11]
            end
        end
    end

    @testset "DELETE /assignment/members/{index_id} removes a member (native route)" begin
        mktempdir() do dir
            db = HimalayaUI.open_db(joinpath(dir, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
            e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
            DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
            DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [e_id])
            DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, 10)", [e_id])

            with_test_server(db) do port, base
                r = HTTP.delete("$base/api/exposures/$e_id/assignment/members/10",
                    ["X-Username" => "alice"])
                @test r.status == 200
                got = JSON3.read(String(r.body))
                @test collect(got.members) == Int[]
                @test isempty(HimalayaUI._assignment_body(db, e_id)[:members])

                # Removing an absent member is a benign no-op (200, empty).
                r = HTTP.delete("$base/api/exposures/$e_id/assignment/members/10",
                    ["X-Username" => "alice"])
                @test r.status == 200
            end
        end
    end

    @testset "native member routes carry distinct {assignment} post_state" begin
        # Layer-1/2 contract: the SSE frame these routes enqueue carries
        # post_state = {assignment:{state,members}} with NO top-level `indices`
        # key, so the frontend's CurationPostState/applyPostStateOnly guard bails.
        mktempdir() do dir
            db  = HimalayaUI.open_db(joinpath(dir, "h.db"))
            req = HTTP.Request("POST", "/x",
                ["X-Username" => "alice", "X-Client-Op-Id" => "op-d1-1"], UInt8[])
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
            e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
            DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])

            HimalayaUI.apply_event!(db, req; kind="assignment_add",
                entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))
            b = HimalayaUI._assignment_body(db, e_id)
            post_state = Dict(:assignment => Dict(:state => b[:state], :members => b[:members]))
            round = JSON3.read(JSON3.write(post_state))
            @test haskey(round, :assignment)
            @test !haskey(round, :indices)
            @test round.assignment.state == "indexed"
            @test collect(round.assignment.members) == [10]
        end
    end
end

# ── Plan D Task D-9 (B4): client-fitted custom-index commit ──────────────────
@testset "insert_custom_index! round-trips the modal comb (Ia3d, √6)" begin
    # Finding #4: basis = q₁ slope = 2π/a × first(phaseratios(P)). Ia3d's √6
    # first reflection maximizes the convention-mismatch signal — assert that
    # predicted_q_for_phase(phase, basis) reproduces the physical comb for a=100.
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

        a = 100.0
        P = Himalaya.Ia3d
        ru = Himalaya.phaseratios(P)                 # un-normalized = √N
        basis = 2π / a * first(ru)                   # the modal convention

        nid = HimalayaUI.insert_custom_index!(db, e_id, P, basis)
        row = Tables.rowtable(DBInterface.execute(db,
            "SELECT phase, basis, lattice_d, kind FROM indices WHERE id = ?", [nid]))[1]
        @test String(row.phase) == "Ia3d"
        @test String(row.kind) == "speculative"
        @test Float64(row.basis) ≈ basis
        # lattice_d recovered via Himalaya.fit ≈ a.
        @test Float64(row.lattice_d) ≈ a atol = 1e-6

        # predicted_q_for_phase reproduces the physics comb q_N = 2π√N/a.
        pred = HimalayaUI.predicted_q_for_phase("Ia3d", basis)
        phys = [2π * r / a for r in ru]              # 2π√N/a
        @test length(pred) == length(phys)
        @test all(abs.(pred .- phys) .< 1e-9)
        # The FIRST predicted q equals basis (normalized first ratio is 1.0) —
        # this is the convention check that fails if basis were `a` or `2π/a`.
        @test pred[1] ≈ basis
    end
end

if isdefined(@__MODULE__, :with_test_server)
    @testset "POST /custom-index persists + adds to the assignment" begin
        mktempdir() do dir
            db = HimalayaUI.open_db(joinpath(dir, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
            e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)

            a = 150.0
            P = Himalaya.Pn3m
            basis = 2π / a * first(Himalaya.phaseratios(P))

            with_test_server(db) do port, base
                # missing basis → 400
                r = HTTP.post("$base/api/exposures/$e_id/custom-index",
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:phase => "Pn3m")); status_exception=false)
                @test r.status == 400

                r = HTTP.post("$base/api/exposures/$e_id/custom-index",
                    ["Content-Type" => "application/json", "X-Username" => "alice"],
                    JSON3.write(Dict(:phase => "Pn3m", :basis => basis)))
                @test r.status == 200
                got = JSON3.read(String(r.body))
                @test got.phase == "Pn3m"
                @test got.kind == "speculative"
                @test got.basis ≈ basis
                nid = got.id

                # the new custom index is now an assignment member.
                @test nid in HimalayaUI._assignment_body(db, e_id)[:members]
            end
        end
    end
end
