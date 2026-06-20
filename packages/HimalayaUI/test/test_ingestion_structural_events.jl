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

    @testset "get_loads_rollup returns nested Load→Sample→Exposure tree (§8.8 shape)" begin
        db = fresh_db()
        (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

        rollup = HimalayaUI.get_loads_rollup(db, exp_id)
        @test length(rollup) == 1
        load = first(rollup)
        # Load-level field names match the §8.8 contract exactly.
        @test load.load_id == load_id
        @test haskey(load, :load_index)
        @test haskey(load, :session_id)
        @test haskey(load, :start_time)
        @test haskey(load, :end_time)
        @test haskey(load, :frame_count)
        @test haskey(load, :note)
        @test length(load.samples) == 2

        sample_ids = Set(s.sample_id for s in load.samples)
        @test s1_id in sample_ids && s2_id in sample_ids
        for s in load.samples
            # Sample-level field names match the §8.8 LoadSample contract exactly.
            @test haskey(s, :name)
            @test haskey(s, :slot_index)
            @test haskey(s, :grouping_source)
            @test haskey(s, :name_source)
            @test haskey(s, :merged_into_id)
            @test haskey(s, :flag)         # GroupingFlag or nothing
            @test length(s.exposures) == 1
            ex = first(s.exposures)
            # Exposure-leaf field names match the §8.8 LoadExposure contract exactly.
            @test haskey(ex, :id)
            @test haskey(ex, :filename)
            @test haskey(ex, :horizontal_position)
            @test haskey(ex, :timestamp)
        end
    end

    @testset "get_loads_rollup attaches a present flag, and suppresses a dismissed one" begin
        # Phase B Task 12's derive_sample_flags is the authority for `flag`.
        # Phase B is a documented prerequisite; derive_sample_flags must exist on-branch.
        #
        # Use a fresh seed (not seed_two_samples) so ALL exposures on s1 have an
        # explicit horizontal_position — the grouping.jl SplitFlag path calls
        # Float64(p) which errors on `missing` (the SQLite NULL representation).
        db = fresh_db()
        exp_id = HimalayaUI.create_experiment!(db; path="/d", data_dir="/d", analysis_dir="/a")
        load_id = DBInterface.lastrowid(DBInterface.execute(db,
            "INSERT INTO loads (experiment_id, load_index) VALUES (?, 1)", [exp_id]))
        s1_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, load_id=load_id,
            slot_index=1, name="HA85 (S01P01)")
        s2_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, load_id=load_id,
            slot_index=2, name="HA85 (S01P02)")
        # s2 gets one exposure (position-less is fine; s2 has only one so no gap check).
        HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s2_id,
            filename="f02.tif")
        # s1 gets two exposures with a > 0.5 mm gap → SplitFlag.
        # Explicit frame_no ensures ORDER BY frame_no, id produces 0.0 mm before 2.0 mm.
        HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s1_id,
            filename="f01a.tif", horizontal_position=0.0, frame_no=10)
        HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s1_id,
            filename="f01b.tif", horizontal_position=2.0, frame_no=11)

        # Before any dismiss: derive_sample_flags emits the split flag for s1.
        r0 = HimalayaUI.get_loads_rollup(db, exp_id)
        s1_before = first(filter(s -> s.sample_id == s1_id, first(r0).samples))
        @test s1_before.flag !== nothing      # proves the flag is really there to suppress

        # Record a non-undone grouping_flag_dismissed event on s1.
        req = user_req()
        SQLite.transaction(db) do
            HimalayaUI.apply_event!(HimalayaUI.InTransaction(), db, req;
                kind = "grouping_flag_dismissed", entity_type = "sample", entity_id = s1_id,
                payload = Dict(:flag_kind => "split", :experiment_id => exp_id))
        end

        # After the dismiss: the same flag is suppressed to nothing.
        r1 = HimalayaUI.get_loads_rollup(db, exp_id)
        s1_after = first(filter(s -> s.sample_id == s1_id, first(r1).samples))
        @test s1_after.flag === nothing       # proves the dismiss is what cleared it
    end

    @testset "get_loads_rollup tolerates NULL horizontal_position (no Float64(missing) throw)" begin
        # Real PRP-partial data: a sample whose exposures mix a measured position
        # with a NULL one. SQLite returns NULL as `missing`, and derive_sample_flags'
        # `=== nothing` gap-guard does NOT catch `missing` — so Float64(missing)
        # would throw and 500 the grouping page. This pins the boundary.
        db = fresh_db()
        (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)
        HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s1_id,
            filename="p1.tif", horizontal_position=1.0, frame_no=10)
        HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s1_id,
            filename="p2.tif", frame_no=11)   # horizontal_position defaults to NULL

        r = HimalayaUI.get_loads_rollup(db, exp_id)   # must not throw
        s1 = first(filter(s -> s.sample_id == s1_id, first(r).samples))
        @test s1.flag === nothing             # a lone positioned/NULL pair yields no split
    end

    @testset "PATCH /api/samples/{id}/name — renames and records event" begin
        db = fresh_db()
        (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

        # Drive the real registered route via in-process dispatch (no socket).
        # with_inproc_routes binds `db` so the handler's current_db() resolves to it,
        # ensures routes are registered, and passes a call closure. Do NOT call
        # register_grouping_routes!() here — registration happens inside
        # with_inproc_routes → _ensure_inproc_routes! → register_routes!.
        with_inproc_routes(db) do call
            resp = call("PATCH", "/api/samples/$(s1_id)/name";
                headers = ["Content-Type"  => "application/json",
                           "X-Username"    => "alice",
                           "X-Client-Op-Id" => "test-op-rename-1"],
                body = Vector{UInt8}(JSON3.write(Dict(:name => "JC C04 (S01P01)"))))
            @test resp.status == 200

            row = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT name, name_source FROM samples WHERE id = ?", [s1_id])))
            @test String(row.name) == "JC C04 (S01P01)"
            @test String(row.name_source) == "user"

            # A sample_renamed event was written.
            evts = Tables.rowtable(DBInterface.execute(db,
                "SELECT action FROM user_actions WHERE entity_type='sample' AND entity_id=?", [s1_id]))
            @test any(r -> String(r.action) == "sample_renamed", evts)
        end
    end

    # Folded from Task 2: pin that rebuild_views_from_log! tolerates entity_type="sample".
    # rebuild_views_from_log! (events.jl:1148) already declares entity_type::String = "exposure"
    # and threads it into the WHERE clause — no events.jl edit required. This assertion
    # prevents a future refactor from silently breaking the sample-partition path.
    @testset "sample_renamed — event durable + rebuild_views_from_log! tolerates entity_type='sample'" begin
        db = fresh_db()
        (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)
        req = user_req()

        result = SQLite.transaction(db) do
            HimalayaUI.apply_event!(HimalayaUI.InTransaction(), db, req;
                kind        = "sample_renamed",
                entity_type = "sample",
                entity_id   = s1_id,
                payload     = Dict(:name => "JC C04 (S01P01)", :experiment_id => exp_id))
        end
        @test result.event_id > 0

        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT action, entity_type, entity_id FROM user_actions WHERE id = ?",
            [result.event_id])))
        @test String(row.action) == "sample_renamed"
        @test String(row.entity_type) == "sample"
        @test Int(row.entity_id) == s1_id

        # rebuild_views_from_log! with entity_type="sample" must not throw.
        @test_nowarn HimalayaUI.rebuild_views_from_log!(db, s1_id; entity_type="sample")
    end

    @testset "POST /api/exposures/{id}/move — moves exposure + records event" begin
        db = fresh_db()
        (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

        with_inproc_routes(db) do call
            resp = call("POST", "/api/exposures/$(e1_id)/move";
                headers = ["Content-Type"  => "application/json",
                           "X-Username"    => "alice",
                           "X-Client-Op-Id" => "test-op-move-1"],
                body = Vector{UInt8}(JSON3.write(Dict(:sample_id => s2_id))))
            @test resp.status == 200

            row = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT sample_id FROM exposures WHERE id = ?", [e1_id])))
            @test Int(row.sample_id) == s2_id

            evts = Tables.rowtable(DBInterface.execute(db,
                "SELECT action FROM user_actions WHERE entity_type='exposure' AND entity_id=?", [e1_id]))
            @test any(r -> String(r.action) == "exposure_moved", evts)
        end
    end

    @testset "POST /api/exposures/{id}/move — rejects cross-experiment move" begin
        db = fresh_db()
        (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)
        # Create a sample in a different experiment.
        exp2_id = HimalayaUI.create_experiment!(db; path="/d2", data_dir="/d2", analysis_dir="/a2")
        s_other = HimalayaUI.create_sample!(db; experiment_id=exp2_id, name="Other")

        with_inproc_routes(db) do call
            resp = call("POST", "/api/exposures/$(e1_id)/move";
                headers = ["Content-Type"  => "application/json",
                           "X-Username"    => "alice",
                           "X-Client-Op-Id" => "test-op-move-2"],
                body = Vector{UInt8}(JSON3.write(Dict(:sample_id => s_other))))
            @test resp.status == 422   # cross-experiment forbidden
        end
    end

    @testset "POST /api/samples/{id}/merge — basic re-point + retire" begin
        db = fresh_db()
        (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

        with_inproc_routes(db) do call
            # Merge s2 (loser) into s1 (survivor).
            resp = call("POST", "/api/samples/$(s2_id)/merge";
                headers = ["Content-Type"  => "application/json",
                           "X-Username"    => "alice",
                           "X-Client-Op-Id" => "test-op-merge-1"],
                body = Vector{UInt8}(JSON3.write(Dict(:survivor_id => s1_id))))
            @test resp.status == 200

            # e2 (was in s2) now belongs to s1.
            row = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT sample_id FROM exposures WHERE id = ?", [e2_id])))
            @test Int(row.sample_id) == s1_id

            # s2 is retired, not deleted.
            s2_row = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT merged_into_id FROM samples WHERE id = ?", [s2_id])))
            @test Int(s2_row.merged_into_id) == s1_id

            # exposure_moved events were recorded for the moved exposures.
            evts = Tables.rowtable(DBInterface.execute(db,
                "SELECT action FROM user_actions WHERE entity_type='exposure' AND entity_id=?", [e2_id]))
            @test any(r -> String(r.action) == "exposure_moved", evts)
        end
    end

    @testset "POST /api/samples/{id}/merge — series_samples dedup: drops loser membership" begin
        db = fresh_db()
        (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

        # Create a series containing both s1 and s2.
        series_id = DBInterface.lastrowid(DBInterface.execute(db,
            "INSERT INTO series DEFAULT VALUES"))
        DBInterface.execute(db,
            "INSERT INTO series_samples (series_id, sample_id, position) VALUES (?, ?, 1)",
            [series_id, s1_id])
        DBInterface.execute(db,
            "INSERT INTO series_samples (series_id, sample_id, position) VALUES (?, ?, 2)",
            [series_id, s2_id])

        with_inproc_routes(db) do call
            resp = call("POST", "/api/samples/$(s2_id)/merge";
                headers = ["Content-Type"  => "application/json",
                           "X-Username"    => "alice",
                           "X-Client-Op-Id" => "test-op-merge-dedup-ss"],
                body = Vector{UInt8}(JSON3.write(Dict(:survivor_id => s1_id))))
            @test resp.status == 200

            # s2's membership in the series where s1 already appears must be dropped.
            surviving = Tables.rowtable(DBInterface.execute(db,
                "SELECT sample_id FROM series_samples WHERE series_id = ?", [series_id]))
            remaining_ids = Set(Int(r.sample_id) for r in surviving)
            @test s2_id ∉ remaining_ids   # loser's duplicate membership dropped
            @test s1_id in remaining_ids  # survivor's membership intact
        end
    end

    @testset "POST /api/samples/{id}/merge — sample_tags dedup: survivor wins" begin
        db = fresh_db()
        (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)

        # Both samples have a tag with key "condition"; values differ.
        DBInterface.execute(db,
            "INSERT INTO sample_tags (sample_id, key, value) VALUES (?, 'condition', 'acid')", [s1_id])
        DBInterface.execute(db,
            "INSERT INTO sample_tags (sample_id, key, value) VALUES (?, 'condition', 'base')", [s2_id])
        # s2 has a unique tag that should survive the merge.
        DBInterface.execute(db,
            "INSERT INTO sample_tags (sample_id, key, value) VALUES (?, 'ph', '8.0')", [s2_id])

        with_inproc_routes(db) do call
            resp = call("POST", "/api/samples/$(s2_id)/merge";
                headers = ["Content-Type"  => "application/json",
                           "X-Username"    => "alice",
                           "X-Client-Op-Id" => "test-op-merge-tags"],
                body = Vector{UInt8}(JSON3.write(Dict(:survivor_id => s1_id))))
            @test resp.status == 200

            tags = Tables.rowtable(DBInterface.execute(db,
                "SELECT key, value FROM sample_tags WHERE sample_id = ?", [s1_id]))
            tag_map = Dict(String(t.key) => String(t.value) for t in tags)
            @test tag_map["condition"] == "acid"  # survivor wins on collision
            @test tag_map["ph"] == "8.0"          # loser's unique tag transferred
        end
    end

    @testset "POST /api/samples/{id}/split — creates new sample + moves exposures" begin
        db = fresh_db()
        (exp_id, load_id, s1_id, s2_id, e1_id, e2_id) = seed_two_samples(db)
        # Add a second exposure to s1 so we have two to split.
        e3_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id,
            sample_id=s1_id, filename="f03.tif")

        with_inproc_routes(db) do call
            # Split: move e3 out of s1 into a new sample.
            resp = call("POST", "/api/samples/$(s1_id)/split";
                headers = ["Content-Type"  => "application/json",
                           "X-Username"    => "alice",
                           "X-Client-Op-Id" => "test-op-split-1"],
                body = Vector{UInt8}(JSON3.write(Dict(:exposure_ids => [e3_id], :name => "HA85 (S01P01b)"))))
            @test resp.status == 201

            resp_body = JSON3.read(String(resp.body))
            new_sample_id = Int(resp_body.new_sample_id)

            # e3 now belongs to the new sample.
            row = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT sample_id FROM exposures WHERE id = ?", [e3_id])))
            @test Int(row.sample_id) == new_sample_id

            # e1 still belongs to s1.
            row2 = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT sample_id FROM exposures WHERE id = ?", [e1_id])))
            @test Int(row2.sample_id) == s1_id

            # New sample has the given name.
            new_row = first(Tables.rowtable(DBInterface.execute(db,
                "SELECT name, experiment_id FROM samples WHERE id = ?", [new_sample_id])))
            @test String(new_row.name) == "HA85 (S01P01b)"
            @test Int(new_row.experiment_id) == exp_id

            # exposure_moved event recorded for each moved exposure.
            evts = Tables.rowtable(DBInterface.execute(db,
                "SELECT action FROM user_actions WHERE entity_type='exposure' AND entity_id=?", [e3_id]))
            @test any(r -> String(r.action) == "exposure_moved", evts)
        end
    end
end
