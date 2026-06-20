# packages/HimalayaUI/test/test_ingestion_structural_events.jl
using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# Standalone-run support: under runtests.jl the in-process harness
# (with_inproc_routes, open_prepared_clone, bind_db!) is already loaded
# (test_http.jl transitively includes test_inproc.jl + test_template_db.jl);
# on a standalone run pull it in. with_test_server is defined by test_http.jl
# so its presence is a reliable sentinel that all helpers are loaded.
# Mirrors test_assignment_reattach.jl:28–30.
if !isdefined(@__MODULE__, :with_test_server)
    include("test_http.jl")
end

# ── helpers ─────────────────────────────────────────────────────────────────

"Open a fresh DB with the full current schema (template-clone; fast).
Uses block-less `mktempdir()` so the backing dir is cleaned at PROCESS exit,
not when this function returns — the returned handle must outlive the call
(the block form `mktempdir() do … end` would rm the dir before any query runs)."
function fresh_db()
    open_prepared_clone(mktempdir())
end

"Seed one experiment + one load + two samples + two exposures (one per sample).
Caller already holds `db`; returns (exp_id, load_id, s1_id, s2_id, e1_id, e2_id)."
function seed_two_samples(db)
    exp_id = HimalayaUI.create_experiment!(db; path="/d", data_dir="/d", analysis_dir="/a")
    # load_id: Phase A adds the loads table; if running in isolation seed it manually.
    load_id = DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO loads (experiment_id, load_index) VALUES (?, 1)", [exp_id]))
    s1_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, load_id=load_id,
        slot_index=1, name="HA85 (S01P01)")
    s2_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, load_id=load_id,
        slot_index=2, name="HA85 (S01P02)")
    e1_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id,
        sample_id=s1_id, filename="f01.tif")
    e2_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id,
        sample_id=s2_id, filename="f02.tif")
    (exp_id, load_id, s1_id, s2_id, e1_id, e2_id)
end

"POST-like request with X-Username, for direct apply_event! unit tests
(Tasks 1–2). Route tests drive in-process HTTP via with_inproc_routes instead."
user_req(name="alice") = HTTP.Request("POST", "/x", ["X-Username" => name], UInt8[])

@testset "structural events (Phase D)" begin
    # task testsets append below

    @testset "exposure_moved — dispatcher writes exposures.sample_id" begin
        db = fresh_db()
        (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)
        req = user_req()

        # Move e1 (currently in s1) to s2. Use the InTransaction variant — the
        # broadcasting default-method variant calls broadcast_event! → current_db()
        # (events.jl:1106) on a non-nothing user_id, which errors with no server
        # bound. InTransaction directly exercises the dispatcher and does NOT
        # broadcast (events.jl:101). It assumes an open tx, so wrap it.
        result = SQLite.transaction(db) do
            HimalayaUI.apply_event!(HimalayaUI.InTransaction(), db, req;
                kind = "exposure_moved",
                entity_type = "exposure",
                entity_id = e1_id,
                payload = Dict(:sample_id => s2_id,
                               :from_sample_id => s1_id,
                               :experiment_id => exp_id))
        end
        @test result.event_id > 0

        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id FROM exposures WHERE id = ?", [e1_id])))
        @test Int(row.sample_id) == s2_id

        # The dispatcher ignores the extra payload keys — it reads only sample_id.
        # Confirm the durable payload round-trips all three fields for the SSE arm.
        pl = JSON3.read(String(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT payload FROM user_actions WHERE id = ?", [result.event_id]))).payload))
        @test Int(pl.sample_id) == s2_id
        @test Int(pl.from_sample_id) == s1_id
        @test Int(pl.experiment_id) == exp_id
    end

    @testset "exposure_moved — rebuild_views_from_log! round-trip" begin
        db = fresh_db()
        (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)
        req = user_req()

        SQLite.transaction(db) do
            HimalayaUI.apply_event!(HimalayaUI.InTransaction(), db, req;
                kind = "exposure_moved", entity_type = "exposure", entity_id = e1_id,
                payload = Dict(:sample_id => s2_id,
                               :from_sample_id => s1_id, :experiment_id => exp_id))
        end

        # Reset view state (undo the dispatcher's write) and replay from log.
        DBInterface.execute(db, "UPDATE exposures SET sample_id = ? WHERE id = ?", [s1_id, e1_id])
        @test Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id FROM exposures WHERE id = ?", [e1_id]))).sample_id) == s1_id

        HimalayaUI.rebuild_views_from_log!(db, e1_id; entity_type="exposure")
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id FROM exposures WHERE id = ?", [e1_id])))
        @test Int(row.sample_id) == s2_id
    end

    @testset "samples.merged_into_id column exists after migration" begin
        db = fresh_db()
        col_names = String.(getproperty.(Tables.rowtable(
            DBInterface.execute(db, "PRAGMA table_info(samples)")), :name))
        @test "merged_into_id" in col_names
    end

    @testset "retire_sample! sets merged_into_id" begin
        db = fresh_db()
        (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

        HimalayaUI.retire_sample!(db, s2_id; merged_into_id=s1_id)

        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT merged_into_id FROM samples WHERE id = ?", [s2_id])))
        @test Int(row.merged_into_id) == s1_id

        # The loser row still exists (no hard delete).
        @test !isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM samples WHERE id = ?", [s2_id])))
    end
end
