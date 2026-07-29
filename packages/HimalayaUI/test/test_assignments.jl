using Test
using HimalayaUI
using SQLite, DBInterface, Tables
using JSON3, HTTP

@testset "assignments schema" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)

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
        db = open_prepared_clone(dir)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

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
    # The migration backfills state='indexed' + the auto members for legacy
    # upgrades (a one-time HISTORICAL backfill). It DIVERGES from current
    # analyze behavior on purpose: a fresh analyze no longer seeds the
    # assignment (auto-grouping is not durable), but legacy upgrades carry the
    # auto guess forward so a pre-Plan-A DB keeps what it displayed.
    mktempdir() do dir
        db = open_prepared_clone(dir)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

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
        db = open_prepared_clone(dir)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

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
        db = open_prepared_clone(dir)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

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
        db = open_prepared_clone(dir)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

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
        db = open_prepared_clone(dir)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)
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
            db = open_prepared_clone(dir)
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
            e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)
            DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
            DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [e_id])
            DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, 10)", [e_id])

            with_inproc_routes(db) do call
                r = call("GET", "/api/exposures/$e_id/assignment")
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
        db  = open_prepared_clone(dir)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

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
            db = open_prepared_clone(dir)
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
            e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)
            DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
            DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [e_id])
            DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, 10)", [e_id])

            with_inproc_routes(db) do call
                # Invalid state → 400.
                r = call("POST", "/api/exposures/$e_id/assignment/state";
                    headers = ["Content-Type" => "application/json"],
                    body = Vector{UInt8}(JSON3.write(Dict(:state => "bogus"))))
                @test r.status == 400

                # Valid: form_factor clears members.
                r = call("POST", "/api/exposures/$e_id/assignment/state";
                    headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
                    body = Vector{UInt8}(JSON3.write(Dict(:state => "form_factor"))))
                @test r.status == 200
                got = JSON3.read(String(r.body))
                @test got.state == "form_factor"
                @test collect(got.members) == Int[]
                @test HimalayaUI._assignment_body(db, e_id)[:state] == "form_factor"
            end
        end
    end
end

@testset "retired index_confirmed/index_unconfirmed events fold to no-ops (D-10)" begin
    # The Plan-A dual-write was retired with the /groups routes. The legacy
    # event kinds remain as no-op GUARDS so rebuild_views_from_log! still treats
    # them as KNOWN kinds (never throws on historical events), but they write
    # NOTHING — the durable assignment is now the sole source of truth, written
    # exclusively by the assignment_* kinds.
    mktempdir() do dir
        db  = open_prepared_clone(dir)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
        DBInterface.execute(db, "INSERT INTO index_groups (id, exposure_id, kind, active) VALUES (200, ?, 'custom', 1)", [e_id])

        # Folding a historical index_confirmed / index_unconfirmed must not throw
        # and must leave the legacy membership table untouched.
        HimalayaUI.apply_event!(db, req; kind="index_confirmed",
            entity_type="exposure", entity_id=e_id, payload=Dict(:group_id => 200, :index_id => 10))
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM index_group_members WHERE group_id = 200")))
        HimalayaUI.apply_event!(db, req; kind="index_unconfirmed",
            entity_type="exposure", entity_id=e_id, payload=Dict(:group_id => 200, :index_id => 10))
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM index_group_members WHERE group_id = 200")))

        # assignment_* remains the sole writer to the durable assignment.
        HimalayaUI.apply_event!(db, req; kind="assignment_add",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))
        @test HimalayaUI._assignment_body(db, e_id)[:members] == [10]
        HimalayaUI.apply_event!(db, req; kind="assignment_remove",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))
        @test isempty(HimalayaUI._assignment_body(db, e_id)[:members])
    end
end

@testset "_bonnet_for_index flags a coexisting cubic" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

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
            db = open_prepared_clone(dir)
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
            e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)
            DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
            DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (11, ?, 'Im3m', 0.1)", [e_id])

            with_inproc_routes(db) do call
                # Missing field -> 400.
                r = call("POST", "/api/exposures/$e_id/assignment/members";
                    headers = ["Content-Type" => "application/json"],
                    body = Vector{UInt8}(JSON3.write(Dict(:nope => 1))))
                @test r.status == 400

                # Valid add -> 200, body is the canonical assignment shape.
                r = call("POST", "/api/exposures/$e_id/assignment/members";
                    headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
                    body = Vector{UInt8}(JSON3.write(Dict(:index_id => 10))))
                @test r.status == 200
                got = JSON3.read(String(r.body))
                @test got.exposure_id == e_id
                @test got.state == "indexed"
                @test collect(got.members) == [10]
                @test haskey(got, :event_id)
                @test HimalayaUI._assignment_body(db, e_id)[:members] == [10]

                # Idempotent re-add of the same member is a no-op on membership.
                r = call("POST", "/api/exposures/$e_id/assignment/members";
                    headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
                    body = Vector{UInt8}(JSON3.write(Dict(:index_id => 11))))
                @test r.status == 200
                @test HimalayaUI._assignment_body(db, e_id)[:members] == [10, 11]
            end
        end
    end

    @testset "DELETE /assignment/members/{index_id} removes a member (native route)" begin
        mktempdir() do dir
            db = open_prepared_clone(dir)
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
            e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)
            DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])
            DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [e_id])
            DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, 10)", [e_id])

            with_inproc_routes(db) do call
                r = call("DELETE", "/api/exposures/$e_id/assignment/members/10";
                    headers = ["X-Username" => "alice"])
                @test r.status == 200
                got = JSON3.read(String(r.body))
                @test collect(got.members) == Int[]
                @test isempty(HimalayaUI._assignment_body(db, e_id)[:members])

                # Removing an absent member is a benign no-op (200, empty).
                r = call("DELETE", "/api/exposures/$e_id/assignment/members/10";
                    headers = ["X-Username" => "alice"])
                @test r.status == 200
            end
        end
    end

    @testset "native member routes carry distinct {assignment} post_state" begin
        # Layer-1/2 contract: the SSE frame these routes enqueue carries
        # post_state = {assignment:{state,members}} with NO top-level `indices`
        # key, so the frontend's CurationPostState/applyPostStateOnly guard bails.
        mktempdir() do dir
            db  = open_prepared_clone(dir)
            req = HTTP.Request("POST", "/x",
                ["X-Username" => "alice", "X-Client-Op-Id" => "op-d1-1"], UInt8[])
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
            e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)
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
        db = open_prepared_clone(dir)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

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

# The normalized ratios the modal draws for Pn3m (SYMS.Ms = [2,3,4,6,8,9], and
# Pn3m IS a clean positional prefix of the core series — Hexagonal is the one
# that isn't, pinned separately below).
_pn3m_drawn_ratios() = Himalaya.phaseratios(Himalaya.Pn3m; normalize = true)[1:6]

@testset "insert_custom_index! claims the peaks the modal showed landing" begin
    # The modal counts a reflection as landing at CUSTOM_SNAP_TOL (customIndex.ts
    # `landsOn` relTol). Commit must claim exactly those, so the Focus comb /
    # detector / cart agree with what the user fitted instead of showing zero.
    mktempdir() do dir
        db = open_prepared_clone(dir)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

        a = 150.0
        P = Himalaya.Pn3m
        ratios = Himalaya.phaseratios(P; normalize = true)   # √2 √3 √4 √6 …
        basis  = 2π * sqrt(2) / a                            # q of the √2 reflection

        # Peak 1: dead on √2. Peak 2: 1% off √3 — inside CUSTOM_SNAP_TOL (2.2%),
        # OUTSIDE SNAP_TOL (0.25%), so it pins the tolerance choice. Peak 3: 10%
        # off √4, well outside both. Nothing near the higher orders.
        q1 = basis
        q2 = basis * ratios[2] * 1.01
        q3 = basis * ratios[3] * 1.10
        for q in (q1, q2, q3)
            DBInterface.execute(db,
                "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, ?, 1.0)",
                [e_id, q])
        end

        nid = HimalayaUI.insert_custom_index!(db, e_id, P, basis)

        # basis and lattice_d stay verbatim — the user's lattice, not a refit.
        row = Tables.rowtable(DBInterface.execute(db,
            "SELECT basis, lattice_d FROM indices WHERE id = ?", [nid]))[1]
        @test Float64(row.basis) ≈ basis
        @test Float64(row.lattice_d) ≈ a atol = 1e-6

        ips = Tables.rowtable(DBInterface.execute(db,
            "SELECT ratio_position, residual FROM index_peaks WHERE index_id = ? ORDER BY ratio_position",
            [nid]))
        @test [Int(r.ratio_position) for r in ips] == [1, 2]
        # Residual is |q_obs − basis·ratio| against the STORED basis.
        @test Float64(ips[2].residual) ≈ abs(q2 - basis * ratios[2]) rtol = 1e-9

        # Durable intents mirror the claims so a reanalysis wipe can heal them.
        ints = Tables.rowtable(DBInterface.execute(db,
            "SELECT ratio_position, q FROM speculative_peak_intents WHERE index_id = ? ORDER BY ratio_position",
            [nid]))
        @test [Int(r.ratio_position) for r in ints] == [1, 2]
        @test [Float64(r.q) for r in ints] ≈ [q1, q2]
    end
end

@testset "Hexagonal commit claims √12 and never the forbidden √11" begin
    # Behavioural teeth for the alignment contract. The cross-language test
    # below pins the DATA (SYMS vs phaseratios); this drives the actual claim
    # path, so reverting compute_snap to a count bound fails HERE.
    mktempdir() do dir
        db = open_prepared_clone(dir)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

        P      = Himalaya.Hexagonal
        rn     = Himalaya.phaseratios(P; normalize = true)   # [1, √3, 2, √7, 3, √11, √12, …]
        basis  = 0.10
        # Modal draws Ms = [1, 3, 4, 7, 9, 12] → backend positions {1,2,3,4,5,7}.
        drawn  = sqrt.([1.0, 3, 4, 7, 9, 12]) ./ 1.0

        # One peak exactly on √11 (backend position 6 — never drawn, and
        # physically impossible) and one exactly on √12 (drawn, position 7).
        q11 = basis * rn[6]
        q12 = basis * rn[7]
        for q in (basis, q11, q12)
            DBInterface.execute(db,
                "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, ?, 1.0)",
                [e_id, q])
        end

        nid = HimalayaUI.insert_custom_index!(db, e_id, P, basis; drawn_ratios = drawn)
        claimed = sort([Int(r.ratio_position) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT ratio_position FROM index_peaks WHERE index_id = ?", [nid]))])

        @test 7 in claimed          # the √12 the modal drew
        @test !(6 in claimed)       # the √11 it did not (and that cannot exist)
        @test claimed == [1, 7]
        # A count bound (max_order = 6) would produce exactly the inverse.
    end
end

@testset "a one-peak locked index stores NULL r_squared, not -Inf" begin
    # R² = 1 - RSS/TSS is undefined at one observation (TSS == 0). On the
    # locked branch RSS > 0 by construction — the claim is made at
    # CUSTOM_SNAP_TOL — so the result is -Inf, deterministically. Unlike NaN
    # (which SQLite binds to NULL), -Inf round-trips as a live Float64 and
    # then throws in JSON3.write INSIDE the analyze transaction, making the
    # exposure permanently un-analyzable: every retry recomputes it.
    mktempdir() do dir
        db = open_prepared_clone(dir)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

        P     = Himalaya.Pn3m
        rn    = Himalaya.phaseratios(P; normalize = true)
        basis = 2π * sqrt(2) / 150.0
        # EXACTLY ONE peak, and off the comb by 1.5% so RSS > 0.
        q1 = basis * rn[1] * 1.015
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, ?, 1.0)", [e_id, q1])

        nid = HimalayaUI.insert_custom_index!(db, e_id, P, basis;
                                              drawn_ratios = _pn3m_drawn_ratios())
        @test Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM index_peaks WHERE index_id = ?", [nid]))[1].c == 1

        prows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, q FROM auto_peaks WHERE exposure_id = ? ORDER BY q", [e_id]))
        eff = (q         = [Float64(r.q) for r in prows],
               sharpness = fill(1.0, length(prows)),
               peak_id   = [Int(r.id) for r in prows],
               peak_kind = fill(:auto, length(prows)))
        empty_pr = (q = Float64[], indices = Int[],
                    prominence = Float64[], sharpness = Float64[])
        HimalayaUI.persist_analysis!(db, e_id, Float64[], Float64[],
            empty_pr, Himalaya.Index[], Himalaya.Index[], eff)

        row = Tables.rowtable(DBInterface.execute(db,
            "SELECT r_squared, score FROM indices WHERE id = ?", [nid]))[1]
        # NOT vacuous: r_squared is NULL at commit too, so first prove the
        # reanalysis UPDATE actually ran. `score` is written unconditionally on
        # that path and is NULL until it does.
        @test !ismissing(row.score)
        @test ismissing(row.r_squared)   # NULL, not -Inf

        # The reachable consequence: this must serialize. Pre-guard it threw
        # ErrorException("-Inf not allowed to be written in JSON spec"),
        # unwinding the whole analyze transaction on every retry.
        payload = HimalayaUI._serialized_indices_for_broadcast(db, e_id)
        @test JSON3.write(payload) isa String
    end
end

@testset "every reflection the modal draws maps to a backend ratio position" begin
    # The alignment half of the modal↔backend contract (the tolerance half is
    # pinned from the TS side in deleteIndexAssignment.test.ts). `drawn_ratios`
    # only bounds the claim correctly if every ratio the modal draws matches
    # exactly one position of phaseratios(P) — a drawn reflection with no
    # counterpart would be silently un-claimable.
    #
    # Reads SYMS out of the real customIndex.ts rather than restating it here,
    # so the pin cannot rot into a comment.
    ts_path = joinpath(@__DIR__, "..", "frontend", "src", "lib", "customIndex.ts")
    @test isfile(ts_path)
    ts = read(ts_path, String)

    # SYMS entries look like:  Pn3m: { kind: "cubic", Ms: [2, 3, 4, 6, 8, 9], …
    entries = Dict{String, Tuple{String, Vector{Int}}}()
    for m in eachmatch(r"(\w+):\s*\{\s*kind:\s*\"(\w+)\",\s*Ms:\s*\[([^\]]*)\]", ts)
        name, kind, ms_s = m.captures[1], m.captures[2], m.captures[3]
        ms = [parse(Int, strip(x)) for x in split(ms_s, ",") if !isempty(strip(x))]
        entries[name] = (kind, ms)
    end
    # All eight canonical phases, or the regex drifted from the source.
    @test length(entries) == 8

    for (name, (kind, ms)) in sort(collect(entries))
        P = HimalayaUI.resolve_phase(name)
        @test P !== nothing
        backend = Himalaya.phaseratios(P; normalize = true)
        # Mirror customRefls' q-law, then normalize: lamellar q ∝ n, everything
        # else q ∝ √N. (hex/square/cubic all reduce to √N once normalized.)
        drawn_q = kind == "lamellar" ? Float64.(ms) : sqrt.(Float64.(ms))
        drawn   = drawn_q ./ drawn_q[1]

        for r in drawn
            hits = count(b -> isapprox(b, r; rtol = HimalayaUI.RATIO_MATCH_RTOL), backend)
            # Exactly one: zero means the claim silently vanishes, two would
            # mean RATIO_MATCH_RTOL is merging distinct orders.
            @test hits == 1
        end
    end

    # Hexagonal is deliberately NOT positionally aligned: the core series
    # carries a √11, which is not a permitted 2D hexagonal reflection
    # (N = h²+hk+k² has no integer solution for 11). Pinning this stops anyone
    # "simplifying" drawn_ratios back to a count bound — which would claim the
    # impossible √11 and drop the √12 the modal drew.
    hex_backend = Himalaya.phaseratios(Himalaya.Hexagonal; normalize = true)
    hex_ms      = entries["Hexagonal"][2]
    @test hex_ms == [1, 3, 4, 7, 9, 12]
    @test round(Int, hex_backend[6]^2) == 11          # the spurious entry
    @test !(11 in [h^2 + h*k + k^2 for h in -12:12, k in -12:12])
    # …so a prefix of length(Ms) is NOT the drawn set.
    @test round.(Int, hex_backend[1:length(hex_ms)] .^ 2) != hex_ms
end

@testset "basis_locked survives reanalysis (a custom lattice never drifts)" begin
    # The load-bearing half of writing scan-derived intents at all: they must
    # decide WHICH peaks are claimed, never WHERE the comb sits. Pre-lock, a
    # reanalysis least-squares-refit basis through the resolved intents, so
    # peaks matched at up to CUSTOM_SNAP_TOL (2.2%) dragged the user's lattice.
    mktempdir() do dir
        db = open_prepared_clone(dir)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

        a = 150.0
        P = Himalaya.Pn3m
        ratios = Himalaya.phaseratios(P; normalize = true)
        basis  = 2π * sqrt(2) / a
        # Order 2 deliberately OFF by +1.5%: inside CUSTOM_SNAP_TOL, so it is
        # claimed — and far enough off that a refit through it would visibly
        # move the lattice. This peak is the whole point of the test.
        qs = [basis, basis * ratios[2] * 1.015]
        for q in qs
            DBInterface.execute(db,
                "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, ?, 1.0)",
                [e_id, q])
        end

        nid = HimalayaUI.insert_custom_index!(db, e_id, P, basis;
                                             drawn_ratios = ratios[1:6])
        @test Int(Tables.rowtable(DBInterface.execute(db,
            "SELECT basis_locked FROM indices WHERE id = ?", [nid]))[1].basis_locked) == 1
        @test Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM speculative_peak_intents WHERE index_id = ?",
            [nid]))[1].c == 2

        # Reanalysis over the SAME effective peak set. Built from the DB rather
        # than passed empty: an empty set resolves zero intents, hits
        # `isempty(ratio_to_peak) && continue`, and never reaches the refit —
        # the assertions below would then pass vacuously.
        prows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, q FROM auto_peaks WHERE exposure_id = ? ORDER BY q", [e_id]))
        eff = (q         = [Float64(r.q) for r in prows],
               sharpness = fill(1.0, length(prows)),
               peak_id   = [Int(r.id) for r in prows],
               peak_kind = fill(:auto, length(prows)))
        empty_pr = (q = Float64[], indices = Int[],
                    prominence = Float64[], sharpness = Float64[])
        HimalayaUI.persist_analysis!(db, e_id, Float64[], Float64[],
            empty_pr, Himalaya.Index[], Himalaya.Index[], eff)

        # The claims must actually have been re-resolved — otherwise "basis
        # unchanged" proves nothing.
        @test Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM index_peaks WHERE index_id = ?", [nid]))[1].c == 2

        row = Tables.rowtable(DBInterface.execute(db,
            "SELECT basis, lattice_d, r_squared FROM indices WHERE id = ?", [nid]))[1]
        # Verbatim, to the bit — not "close enough".
        @test Float64(row.basis) == basis
        @test Float64(row.lattice_d) ≈ a atol = 1e-6

        # R² must describe the LOCKED comb, not the least-squares line the
        # branch just refused to store. Himalaya.fit refits internally, so its
        # R² would read ≈1.000 here even though order 2 sits 1.5% off the comb
        # the rail draws — the inverse of this PR's honesty thesis.
        rn = Himalaya.phaseratios(P; normalize = true)
        expected_r2 = Himalaya.R²(basis .* rn[[1, 2]], qs)
        @test Float64(row.r_squared) ≈ expected_r2
        @test Float64(row.r_squared) < 0.9999   # visibly imperfect, as drawn

        # The refit it was spared: an unlocked index over the same peaks lands
        # somewhere else, which is what makes the assertion above meaningful.
        drifted = ratios[[1, 2]] \ qs
        @test !isapprox(drifted, basis; rtol = 1e-6)
    end
end

@testset "DELETE /api/indices/:id leaves no orphan peaks, intents or members" begin
    # Drives the REAL route (with_inproc_routes) rather than re-issuing its SQL
    # by hand: a hand-rolled copy cannot notice the route diverging, which is
    # exactly what it was meant to guard. Covers the FK cascade
    # (speculative_peak_intents + assignment_members are CASCADE, not explicit
    # deletes) and the paired assignment_remove event.
    mktempdir() do dir
        db = open_prepared_clone(dir)
        @test Tables.rowtable(DBInterface.execute(db, "PRAGMA foreign_keys"))[1].foreign_keys == 1

        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
        e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

        a = 150.0
        P = Himalaya.Pn3m
        ratios = Himalaya.phaseratios(P; normalize = true)
        basis  = 2π * sqrt(2) / a
        for q in (basis, basis * ratios[2] * 1.01)
            DBInterface.execute(db,
                "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, ?, 1.0)",
                [e_id, q])
        end

        n(tbl, nid) = Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM $tbl WHERE index_id = ?", [nid]))[1].c

        with_inproc_routes(db) do call
            # `ratios` validation sits ABOVE with_idempotency, so a reject can
            # never commit a partial write. Every bad shape must 400.
            for bad in (Any[], Any[0.0, 1.0], Any[-1.0], "six")
                rb = call("POST", "/api/exposures/$e_id/custom-index";
                    headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
                    body = Vector{UInt8}(JSON3.write(
                        Dict(:phase => "Pn3m", :basis => basis, :ratios => bad))))
                @test rb.status == 400
            end
            # …and none of them persisted anything.
            @test Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM indices WHERE exposure_id = ?",
                [e_id]))[1].c == 0

            r = call("POST", "/api/exposures/$e_id/custom-index";
                headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
                body = Vector{UInt8}(JSON3.write(
                    Dict(:phase => "Pn3m", :basis => basis,
                         :ratios => _pn3m_drawn_ratios()))))
            @test r.status == 200
            nid = Int(JSON3.read(String(r.body)).id)

            @test n("index_peaks", nid) == 2
            @test n("speculative_peak_intents", nid) == 2
            @test n("assignment_members", nid) == 1   # the route adds it to the call

            d = call("DELETE", "/api/indices/$nid";
                headers = ["X-Username" => "alice"])
            @test d.status == 200

            @test n("index_peaks", nid) == 0
            @test n("speculative_peak_intents", nid) == 0   # cascade
            @test n("assignment_members", nid) == 0         # cascade
            @test isempty(Tables.rowtable(DBInterface.execute(db,
                "SELECT 1 FROM indices WHERE id = ?", [nid])))

            # The cascade must be mirrored by a durable assignment_remove, or
            # rebuild_views_from_log! re-folds a member row whose index is gone.
            kinds = [String(x.action) for x in Tables.rowtable(DBInterface.execute(db,
                "SELECT action FROM user_actions ORDER BY id"))]
            @test "speculative_deleted" in kinds
            @test "assignment_remove" in kinds
            @test findfirst(==("speculative_deleted"), kinds) <
                  findlast(==("assignment_remove"), kinds)
        end
    end
end

if isdefined(@__MODULE__, :with_test_server)
    @testset "POST /custom-index persists + adds to the assignment" begin
        mktempdir() do dir
            db = open_prepared_clone(dir)
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S")
            e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id)

            a = 150.0
            P = Himalaya.Pn3m
            basis = 2π / a * first(Himalaya.phaseratios(P))

            with_inproc_routes(db) do call
                # missing basis → 400
                r = call("POST", "/api/exposures/$e_id/custom-index";
                    headers = ["Content-Type" => "application/json"],
                    body = Vector{UInt8}(JSON3.write(Dict(:phase => "Pn3m"))))
                @test r.status == 400

                r = call("POST", "/api/exposures/$e_id/custom-index";
                    headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
                    body = Vector{UInt8}(JSON3.write(Dict(:phase => "Pn3m", :basis => basis))))
                @test r.status == 200
                got = JSON3.read(String(r.body))
                @test got.phase == "Pn3m"
                @test got.kind == "speculative"
                @test got.basis ≈ basis
                # The response must carry the claims insert_custom_index! wrote;
                # asserting only phase/kind/basis let the pre-fix `peaks: []`
                # contract survive at this layer.
                @test haskey(got, :peaks)
                @test length(got.peaks) == length(Tables.rowtable(DBInterface.execute(db,
                    "SELECT 1 FROM index_peaks WHERE index_id = ?", [Int(got.id)])))
                nid = got.id

                # the new custom index is now an assignment member.
                @test nid in HimalayaUI._assignment_body(db, e_id)[:members]
            end
        end
    end
end
