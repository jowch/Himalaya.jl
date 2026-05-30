using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

@testset "apply_event! writes to user_actions and returns event_id" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x",
            ["X-Username" => "alice"], UInt8[])

        result = HimalayaUI.apply_event!(db, req;
            kind = "test_kind", entity_type = "exposure", entity_id = 42,
            payload = Dict(:foo => "bar"))
        eid = result.event_id
        @test eid > 0

        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM user_actions WHERE id = ?", [eid])))
        @test String(row.action) == "test_kind"
        @test String(row.entity_type) == "exposure"
        @test Int(row.entity_id) == 42
        @test JSON3.read(row.payload).foo == "bar"
        # user_id should resolve to alice's id.
        @test !ismissing(row.user_id)
    end
end

@testset "apply_event! with no payload writes NULL payload and skips dispatcher" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "bob"], UInt8[])
        result = HimalayaUI.apply_event!(db, req;
            kind = "no_payload_event", entity_type = "exposure", entity_id = 1)
        eid = result.event_id
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT payload FROM user_actions WHERE id = ?", [eid])))
        @test ismissing(row.payload)
        @test result.view_row_id === nothing
    end
end

@testset "apply_event! with no X-Username writes user_id = NULL" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", Pair{String,String}[], UInt8[])
        result = HimalayaUI.apply_event!(db, req;
            kind = "anon_event", entity_type = "exposure", entity_id = 1,
            payload = Dict(:k => "v"))
        eid = result.event_id
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT user_id FROM user_actions WHERE id = ?", [eid])))
        @test ismissing(row.user_id)
    end
end

@testset "log_action! still works (legacy wrapper)" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        # Existing call sites pass these args; should write a row with the
        # legacy column shape and a payload of {:note => "..."} when given.
        HimalayaUI.log_action!(db, req;
            action = "set_status", entity_type = "exposure", entity_id = 7,
            note = "rejected")
        row = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT action, entity_type, entity_id, payload FROM user_actions
             ORDER BY id DESC LIMIT 1")))
        @test String(row.action) == "set_status"
        @test JSON3.read(row.payload).note == "rejected"
    end
end

@testset "idx_events_by_exposure index exists" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='index' AND name='idx_events_by_exposure'"))
        @test !isempty(rows)
    end
end

@testset "user_actions has payload and undoes_event_id columns after migration" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        cols = [String(c.name) for c in Tables.rowtable(DBInterface.execute(db,
            "PRAGMA table_info(user_actions)"))]
        @test "payload" in cols
        @test "undoes_event_id" in cols
    end
end

@testset "rebuild_views_from_log! reproduces incremental view state" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])

        # Set up an exposure with one auto peak (so peak_excluded has a target).
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, 0.10, 1.0)", [e_id])

        # Apply a sequence of view-producing events.
        r1 = HimalayaUI.apply_event!(db, req;
            kind = "peak_added", entity_type = "exposure", entity_id = e_id,
            payload = Dict(:q => 0.20))
        @test r1.event_id > 0
        @test r1.view_row_id isa Int
        r2 = HimalayaUI.apply_event!(db, req;
            kind = "peak_excluded", entity_type = "exposure", entity_id = e_id,
            payload = Dict(:q => 0.10, :auto_peak_id => 1))
        @test r2.event_id > r1.event_id
        @test r2.view_row_id isa Int

        # Snapshot view state (sort by q for determinism).
        curations_before = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, kind, q FROM peak_curations WHERE exposure_id = ? ORDER BY q",
            [e_id]))

        # Wipe the views (simulating disaster recovery / migration).
        DBInterface.execute(db, "DELETE FROM peak_curations WHERE exposure_id = ?", [e_id])

        # Sanity: wiped.
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM peak_curations WHERE exposure_id = ?", [e_id])))

        # Rebuild from the log.
        HimalayaUI.rebuild_views_from_log!(db, e_id)

        curations_after = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, kind, q FROM peak_curations WHERE exposure_id = ? ORDER BY q",
            [e_id]))

        # Assert view state matches: same number of rows and same (kind, q) tuples.
        @test length(curations_before) == length(curations_after) == 2
        @test Set([(String(c.kind), Float64(c.q)) for c in curations_before]) ==
              Set([(String(c.kind), Float64(c.q)) for c in curations_after])
    end
end

@testset "apply_event! persists client_id when X-Client-Id header present" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        eid    = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="e1")
        req = HTTP.Request("POST", "/x",
            ["X-Username" => "alice", "X-Client-Id" => "tab-xyz"], UInt8[])
        HimalayaUI.apply_event!(db, req;
            kind="noop_test", entity_type="exposure", entity_id=eid, payload=Dict(:k=>1))
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT client_id FROM user_actions ORDER BY id DESC LIMIT 1"))
        @test rows[1].client_id == "tab-xyz"
    end
end

@testset "apply_event! writes NULL client_id when header absent" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        eid    = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="e1")
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        HimalayaUI.apply_event!(db, req;
            kind="noop_test", entity_type="exposure", entity_id=eid, payload=Dict(:k=>1))
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT client_id FROM user_actions ORDER BY id DESC LIMIT 1"))
        @test ismissing(rows[1].client_id)
    end
end

@testset "apply_event! persists client_op_id when X-Client-Op-Id header present" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        eid    = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="e1")
        req = HTTP.Request("POST", "/x",
            ["X-Username" => "alice", "X-Client-Op-Id" => "uuid-abc"], UInt8[])
        result = HimalayaUI.apply_event!(db, req;
            kind="noop_test", entity_type="exposure", entity_id=eid, payload=Dict(:q=>1.0))
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT client_op_id FROM user_actions WHERE id = ?", [result.event_id]))
        @test String(rows[1].client_op_id) == "uuid-abc"
    end
end

@testset "apply_event! writes NULL client_op_id when header absent" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        eid    = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="e1")
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        result = HimalayaUI.apply_event!(db, req;
            kind="noop_test", entity_type="exposure", entity_id=eid, payload=Dict(:q=>1.0))
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT client_op_id FROM user_actions WHERE id = ?", [result.event_id]))
        @test ismissing(rows[1].client_op_id)
    end
end

@testset "rebuild_views_from_log! tolerates rows with NULL/non-NULL client_op_id mix" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="e1")
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, 0.10, 1.0)", [e_id])

        # Mix: one event without X-Client-Op-Id (NULL), one with.
        req_no  = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        req_yes = HTTP.Request("POST", "/x",
            ["X-Username" => "alice", "X-Client-Op-Id" => "uuid-1"], UInt8[])
        HimalayaUI.apply_event!(db, req_no;
            kind="peak_added", entity_type="exposure", entity_id=e_id,
            payload=Dict(:q => 0.20))
        HimalayaUI.apply_event!(db, req_yes;
            kind="peak_added", entity_type="exposure", entity_id=e_id,
            payload=Dict(:q => 0.30))

        # Snapshot live state, wipe, rebuild.
        before = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, kind, q FROM peak_curations WHERE exposure_id = ? ORDER BY q",
            [e_id]))
        DBInterface.execute(db, "DELETE FROM peak_curations WHERE exposure_id = ?", [e_id])
        HimalayaUI.rebuild_views_from_log!(db, Int(e_id))
        after = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, kind, q FROM peak_curations WHERE exposure_id = ? ORDER BY q",
            [e_id]))

        @test length(before) == length(after) == 2
        @test Set([(String(r.kind), Float64(r.q)) for r in before]) ==
              Set([(String(r.kind), Float64(r.q)) for r in after])
    end
end

@testset "rebuild_views_from_log! tolerates rows with NULL client_id" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        eid    = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="e1")
        DBInterface.execute(db, """
            INSERT INTO user_actions (user_id, action, entity_type, entity_id, payload, client_id)
            VALUES (NULL, 'noop_test', 'exposure', ?, '{}', NULL)
        """, [eid])
        @test_nowarn HimalayaUI.rebuild_views_from_log!(db, Int(eid))
    end
end

@testset "migrate_r4_rebase_entity_type! rewrites legacy peak/index entity_type rows" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])

        # Set up minimal data.
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db,
            "INSERT INTO auto_peaks (exposure_id, q, sharpness) VALUES (?, 0.10, 1.0)", [e_id])
        ap_id = Int(DBInterface.lastrowid(DBInterface.execute(db, "SELECT last_insert_rowid()")))

        # Simulate a legacy user_actions row with entity_type='peak'.
        DBInterface.execute(db,
            """INSERT INTO user_actions (action, entity_type, entity_id)
               VALUES ('add_peak', 'peak', ?)""", [ap_id])

        # Run the migration.
        HimalayaUI.migrate_r4_rebase_entity_type!(db)

        # The row should now have entity_type='exposure'.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT entity_type, entity_id FROM user_actions WHERE action = 'add_peak'"))
        @test !isempty(rows)
        @test String(rows[1].entity_type) == "exposure"
        @test Int(rows[1].entity_id) == e_id

        # Idempotent: running again does nothing.
        HimalayaUI.migrate_r4_rebase_entity_type!(db)
        rows2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT entity_type FROM user_actions WHERE action = 'add_peak'"))
        @test String(rows2[1].entity_type) == "exposure"
    end
end

@testset "apply_event! with defer_broadcast=true does NOT fire broadcast" begin
    sub = (id = "t-defer", pending = Channel{String}(8))
    lock(HimalayaUI.SSE_LOCK) do
        push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
    end
    try
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
            HimalayaUI.bind_db!(db)
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
            e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="e1")
            req = HTTP.Request("POST", "/", Pair{String,String}[], UInt8[])
            result = HimalayaUI.apply_event!(db, req;
                kind="noop_test", entity_type="exposure", entity_id=e_id,
                payload=Dict(:q => 1.0),
                defer_broadcast=true)
            @test Base.n_avail(sub.pending) == 0
            @test result.event_id > 0
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM user_actions WHERE id = ?", [result.event_id]))
            @test length(rows) == 1
            SQLite.close(db)
        end
    finally
        lock(HimalayaUI.SSE_LOCK) do
            filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
        end
        close(sub.pending)
        HimalayaUI.SSE_SUBSCRIBERS[] = []
    end
end

@testset "apply_event! defer_broadcast=false (default) preserves existing behavior" begin
    sub = (id = "t-default", pending = Channel{String}(8))
    lock(HimalayaUI.SSE_LOCK) do
        push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
    end
    try
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
            HimalayaUI.bind_db!(db)
            exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
            e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="e1")
            req = HTTP.Request("POST", "/", Pair{String,String}[], UInt8[])
            HimalayaUI.apply_event!(db, req;
                kind="noop_test", entity_type="exposure", entity_id=e_id,
                payload=Dict(:q => 1.0))
            @test Base.n_avail(sub.pending) == 1
            SQLite.close(db)
        end
    finally
        lock(HimalayaUI.SSE_LOCK) do
            filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
        end
        close(sub.pending)
        HimalayaUI.SSE_SUBSCRIBERS[] = []
    end
end

@testset "I2: apply_event! with same (client_op_id, action, entity) returns existing event_id" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        # Seed FK chain.
        DBInterface.execute(db,
            "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e', '/p', '/d', '/a')")
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, name) VALUES (1, 'A1')")
        res = DBInterface.execute(db,
            "INSERT INTO exposures (sample_id, filename) VALUES (1, 'f')")
        exp_id = Int(DBInterface.lastrowid(res))

        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "op-dup"], UInt8[])
        r1 = HimalayaUI.apply_event!(db, req;
            kind="peak_added", entity_type="exposure", entity_id=exp_id,
            payload=Dict(:q => 1.0))
        r2 = HimalayaUI.apply_event!(db, req;
            kind="peak_added", entity_type="exposure", entity_id=exp_id,
            payload=Dict(:q => 1.0))
        @test r1.event_id == r2.event_id

        # Only one user_actions row exists for this op tuple.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM user_actions WHERE client_op_id = 'op-dup'"))
        @test length(rows) == 1

        # Only one peak_curations row exists (the dispatcher's first INSERT).
        crows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM peak_curations WHERE exposure_id = ?", [exp_id]))
        @test length(crows) == 1
        SQLite.close(db)
    end
end

@testset "I2: apply_event!(::InTransaction, ...) participates in caller's transaction" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        DBInterface.execute(db,
            "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e', '/p', '/d', '/a')")
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, name) VALUES (1, 'A1')")
        res = DBInterface.execute(db,
            "INSERT INTO exposures (sample_id, filename) VALUES (1, 'f')")
        exp_id = Int(DBInterface.lastrowid(res))

        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "op-intx"], UInt8[])
        SQLite.transaction(db) do
            r = HimalayaUI.apply_event!(HimalayaUI.InTransaction(), db, req;
                kind="peak_added", entity_type="exposure", entity_id=exp_id,
                payload=Dict(:q => 2.0))
            @test r.event_id > 0
        end
        # Row durable after outer commit.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM user_actions WHERE client_op_id = 'op-intx'"))
        @test length(rows) == 1
        SQLite.close(db)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Compare page (Plan §Phase 1, Task 1.2 + 1.3): comparison_created /
# comparison_submitted / comparison_deleted dispatcher branches.
# ─────────────────────────────────────────────────────────────────────────────

# Helper: build a single member dict with the minimum required fields.
_compare_member_payload(; id=nothing, exposure_id=nothing, display_order::Int=0,
                          band_height::Float64=1.0, y_offset::Float64=0.0,
                          normalization::String="none",
                          color_override=nothing, label_override=nothing,
                          q_window_min=nothing, q_window_max=nothing,
                          peak_display=nothing,
                          snapshot=Dict(:effective_peaks => [],
                                        :confirmed_index => nothing,
                                        :analysis_inputs_hash => "sha256:zero")) =
    Dict{Symbol,Any}(
        :id             => id,
        :exposure_id    => exposure_id,
        :display_order  => display_order,
        :band_height    => band_height,
        :y_offset       => y_offset,
        :normalization  => normalization,
        :color_override => color_override,
        :label_override => label_override,
        :q_window_min   => q_window_min,
        :q_window_max   => q_window_max,
        :peak_display   => peak_display,
        :snapshot       => snapshot,
    )

# Helper: pre-mint a comparisons row at a known id so the dispatcher's
# INSERT/UPDATE-on-existing path runs against an existing row (matches the
# real route handler's two-step pattern).
# Uses NULL placeholders to mirror the post-#67 route — the dispatcher's
# `COALESCE(col, ?)` then stamps real values on first fold.
function _premint_comparison!(db, id::Int)
    DBInterface.execute(db,
        """INSERT INTO comparisons (id, title, content_hash, created_at, updated_at)
           VALUES (?, NULL, NULL, NULL, NULL)""", [id])
    nothing
end

@testset "comparison_created: 0 members" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        HimalayaUI.bind_db!(db)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        _premint_comparison!(db, 1)
        HimalayaUI.apply_event!(db, req;
            kind="comparison_created", entity_type="comparison", entity_id=1,
            payload=Dict(:title => "T1", :description => "d",
                         :forked_from_id => nothing, :forked_at_hash => nothing,
                         :members => []))
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT title, description, content_hash FROM comparisons WHERE id = 1"))
        @test length(rows) == 1
        @test String(rows[1].title) == "T1"
        @test String(rows[1].description) == "d"
        @test startswith(String(rows[1].content_hash), "sha256:")
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparison_members WHERE comparison_id = 1")))
    end
end

@testset "comparison_created: 3 members, content_hash populated, member fields present" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        HimalayaUI.bind_db!(db)
        # Seed exposures so exposure_id FKs are valid.
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e1 = HimalayaUI.create_exposure!(db; sample_id=s_id)
        e2 = HimalayaUI.create_exposure!(db; sample_id=s_id)
        e3 = HimalayaUI.create_exposure!(db; sample_id=s_id)

        _premint_comparison!(db, 1)
        members = [
            _compare_member_payload(exposure_id=e1, display_order=0, band_height=1.0),
            _compare_member_payload(exposure_id=e2, display_order=1, band_height=2.0),
            _compare_member_payload(exposure_id=e3, display_order=2, band_height=3.0),
        ]
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        HimalayaUI.apply_event!(db, req;
            kind="comparison_created", entity_type="comparison", entity_id=1,
            payload=Dict(:title => "Three", :description => nothing,
                         :forked_from_id => nothing, :forked_at_hash => nothing,
                         :members => members))
        rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT exposure_id, display_order, band_height, snapshot, created_by
               FROM comparison_members WHERE comparison_id = 1
               ORDER BY display_order"""))
        @test length(rows) == 3
        @test [Int(r.exposure_id) for r in rows] == [e1, e2, e3]
        @test [Int(r.display_order) for r in rows] == [0, 1, 2]
        @test [Float64(r.band_height) for r in rows] == [1.0, 2.0, 3.0]
        # snapshot column populated (not NULL).
        @test all(!ismissing(r.snapshot) for r in rows)
        # created_by populated from the X-Username actor.
        @test all(!ismissing(r.created_by) for r in rows)

        ch = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT content_hash FROM comparisons WHERE id = 1"))).content_hash
        @test startswith(String(ch), "sha256:")
    end
end

@testset "comparison_submitted: no-op (same state) preserves rows but recomputes hash" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id = HimalayaUI.create_exposure!(db; sample_id=s_id)
        _premint_comparison!(db, 7)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])

        HimalayaUI.apply_event!(db, req;
            kind="comparison_created", entity_type="comparison", entity_id=7,
            payload=Dict(:title => "T", :description => nothing,
                         :forked_from_id => nothing, :forked_at_hash => nothing,
                         :members => [_compare_member_payload(exposure_id=e_id)]))
        member_id = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparison_members WHERE comparison_id = 7"))).id)

        # Submit with the same payload structure (carrying the existing member id).
        HimalayaUI.apply_event!(db, req;
            kind="comparison_submitted", entity_type="comparison", entity_id=7,
            payload=Dict(:title => "T", :description => nothing,
                         :members => [_compare_member_payload(id=member_id, exposure_id=e_id)]))
        # Same single member, same exposure_id.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, exposure_id FROM comparison_members WHERE comparison_id = 7"))
        @test length(rows) == 1
        @test Int(rows[1].id) == member_id
        @test Int(rows[1].exposure_id) == e_id
    end
end

@testset "comparison_submitted: member added (id === nothing)" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e1 = HimalayaUI.create_exposure!(db; sample_id=s_id)
        e2 = HimalayaUI.create_exposure!(db; sample_id=s_id)
        _premint_comparison!(db, 5)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])

        HimalayaUI.apply_event!(db, req;
            kind="comparison_created", entity_type="comparison", entity_id=5,
            payload=Dict(:title => "T", :description => nothing,
                         :forked_from_id => nothing, :forked_at_hash => nothing,
                         :members => [_compare_member_payload(exposure_id=e1, display_order=0)]))
        m1_id = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparison_members WHERE comparison_id = 5"))).id)

        HimalayaUI.apply_event!(db, req;
            kind="comparison_submitted", entity_type="comparison", entity_id=5,
            payload=Dict(:title => "T", :description => nothing,
                         :members => [
                            _compare_member_payload(id=m1_id, exposure_id=e1, display_order=0),
                            _compare_member_payload(id=nothing, exposure_id=e2, display_order=1),
                         ]))
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, display_order FROM comparison_members WHERE comparison_id = 5 ORDER BY display_order"))
        @test length(rows) == 2
        @test [Int(r.exposure_id) for r in rows] == [e1, e2]
    end
end

@testset "comparison_submitted: member removed (DB id absent from payload)" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e1 = HimalayaUI.create_exposure!(db; sample_id=s_id)
        e2 = HimalayaUI.create_exposure!(db; sample_id=s_id)
        _premint_comparison!(db, 8)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])

        HimalayaUI.apply_event!(db, req;
            kind="comparison_created", entity_type="comparison", entity_id=8,
            payload=Dict(:title => "T", :description => nothing,
                         :forked_from_id => nothing, :forked_at_hash => nothing,
                         :members => [
                            _compare_member_payload(exposure_id=e1, display_order=0),
                            _compare_member_payload(exposure_id=e2, display_order=1),
                         ]))
        ids = [Int(r.id) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparison_members WHERE comparison_id = 8 ORDER BY display_order"))]

        HimalayaUI.apply_event!(db, req;
            kind="comparison_submitted", entity_type="comparison", entity_id=8,
            payload=Dict(:title => "T", :description => nothing,
                         :members => [
                            _compare_member_payload(id=ids[1], exposure_id=e1, display_order=0),
                         ]))
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, exposure_id FROM comparison_members WHERE comparison_id = 8"))
        @test length(rows) == 1
        @test Int(rows[1].id) == ids[1]
    end
end

@testset "comparison_submitted: UPDATE writes snapshot unconditionally on same exposure" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e1 = HimalayaUI.create_exposure!(db; sample_id=s_id)
        _premint_comparison!(db, 12)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])

        first_snap = Dict(:effective_peaks => [Dict(:id => 1, :q => 0.1, :intensity => 1.0,
                                                    :sharpness => 0.5, :source => "auto")],
                          :confirmed_index => nothing,
                          :analysis_inputs_hash => "sha256:firstsnap")
        HimalayaUI.apply_event!(db, req;
            kind="comparison_created", entity_type="comparison", entity_id=12,
            payload=Dict(:title => "T", :description => nothing,
                         :forked_from_id => nothing, :forked_at_hash => nothing,
                         :members => [_compare_member_payload(exposure_id=e1, snapshot=first_snap)]))
        m_id = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparison_members WHERE comparison_id = 12"))).id)

        # Same exposure, NEW snapshot — must update unconditionally.
        new_snap = Dict(:effective_peaks => [],
                        :confirmed_index => Dict(:id => 99, :phase => "Pn3m",
                                                 :lattice_d => 12.0, :r_squared => 0.99,
                                                 :ngc => 0.7, :peak_ids => [42]),
                        :analysis_inputs_hash => "sha256:newsnap")
        HimalayaUI.apply_event!(db, req;
            kind="comparison_submitted", entity_type="comparison", entity_id=12,
            payload=Dict(:title => "T", :description => nothing,
                         :members => [_compare_member_payload(id=m_id, exposure_id=e1,
                                                               snapshot=new_snap)]))
        snap_str = String(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT snapshot FROM comparison_members WHERE id = ?", [m_id]))).snapshot)
        snap_obj = JSON3.read(snap_str)
        @test String(snap_obj.analysis_inputs_hash) == "sha256:newsnap"
        @test snap_obj.confirmed_index !== nothing
        @test String(snap_obj.confirmed_index.phase) == "Pn3m"
    end
end

@testset "comparison_submitted: reorder updates display_order" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e1 = HimalayaUI.create_exposure!(db; sample_id=s_id)
        e2 = HimalayaUI.create_exposure!(db; sample_id=s_id)
        _premint_comparison!(db, 3)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])

        HimalayaUI.apply_event!(db, req;
            kind="comparison_created", entity_type="comparison", entity_id=3,
            payload=Dict(:title => "T", :description => nothing,
                         :forked_from_id => nothing, :forked_at_hash => nothing,
                         :members => [
                            _compare_member_payload(exposure_id=e1, display_order=0),
                            _compare_member_payload(exposure_id=e2, display_order=1),
                         ]))
        rows0 = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, exposure_id, display_order FROM comparison_members WHERE comparison_id = 3 ORDER BY display_order"))
        m1, m2 = Int(rows0[1].id), Int(rows0[2].id)

        # Swap order.
        HimalayaUI.apply_event!(db, req;
            kind="comparison_submitted", entity_type="comparison", entity_id=3,
            payload=Dict(:title => "T", :description => nothing,
                         :members => [
                            _compare_member_payload(id=m1, exposure_id=e1, display_order=1),
                            _compare_member_payload(id=m2, exposure_id=e2, display_order=0),
                         ]))
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, display_order FROM comparison_members WHERE comparison_id = 3 ORDER BY id"))
        @test Int(rows[1].id) == m1
        @test Int(rows[1].display_order) == 1
        @test Int(rows[2].id) == m2
        @test Int(rows[2].display_order) == 0
    end
end

@testset "comparison_submitted: zero-row UPDATE errors (member id not for this comparison)" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        HimalayaUI.bind_db!(db)
        _premint_comparison!(db, 4)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        # Try to UPDATE a member that doesn't exist.
        @test_throws Exception HimalayaUI.apply_event!(db, req;
            kind="comparison_submitted", entity_type="comparison", entity_id=4,
            payload=Dict(:title => "T", :description => nothing,
                         :members => [_compare_member_payload(id=999_999, exposure_id=nothing)]))
    end
end

@testset "comparison_deleted cascades to members and messages" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e1 = HimalayaUI.create_exposure!(db; sample_id=s_id)
        _premint_comparison!(db, 11)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        HimalayaUI.apply_event!(db, req;
            kind="comparison_created", entity_type="comparison", entity_id=11,
            payload=Dict(:title => "T", :description => nothing,
                         :forked_from_id => nothing, :forked_at_hash => nothing,
                         :members => [_compare_member_payload(exposure_id=e1)]))
        DBInterface.execute(db,
            "INSERT INTO comparison_messages (comparison_id, body) VALUES (?, ?)",
            [11, "hi"])
        HimalayaUI.apply_event!(db, req;
            kind="comparison_deleted", entity_type="comparison", entity_id=11,
            payload=Dict(:reason => "test"))
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparisons WHERE id = 11")))
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparison_members WHERE comparison_id = 11")))
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparison_messages WHERE comparison_id = 11")))
    end
end

@testset "rebuild_views_from_log!(entity_type=\"comparison\") reproduces live state" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        HimalayaUI.bind_db!(db)
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e1 = HimalayaUI.create_exposure!(db; sample_id=s_id)
        e2 = HimalayaUI.create_exposure!(db; sample_id=s_id)
        e3 = HimalayaUI.create_exposure!(db; sample_id=s_id)
        cmp_id = 21
        _premint_comparison!(db, cmp_id)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])

        HimalayaUI.apply_event!(db, req;
            kind="comparison_created", entity_type="comparison", entity_id=cmp_id,
            payload=Dict(:title => "T", :description => "d",
                         :forked_from_id => nothing, :forked_at_hash => nothing,
                         :members => [
                            _compare_member_payload(exposure_id=e1, display_order=0),
                            _compare_member_payload(exposure_id=e2, display_order=1),
                         ]))
        ids = [Int(r.id) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM comparison_members WHERE comparison_id = ? ORDER BY display_order",
            [cmp_id]))]
        HimalayaUI.apply_event!(db, req;
            kind="comparison_submitted", entity_type="comparison", entity_id=cmp_id,
            payload=Dict(:title => "T2", :description => "d2",
                         :members => [
                            _compare_member_payload(id=ids[1], exposure_id=e1, display_order=0),
                            _compare_member_payload(id=nothing, exposure_id=e3, display_order=1),
                         ]))

        # Snapshot live state.
        live_cmp = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT title, description, content_hash FROM comparisons WHERE id = ?", [cmp_id])))
        live_members = Tables.rowtable(DBInterface.execute(db,
            """SELECT id, exposure_id, display_order, snapshot
               FROM comparison_members WHERE comparison_id = ? ORDER BY id""", [cmp_id]))

        # Wipe view tables for this comparison and replay.
        DBInterface.execute(db, "DELETE FROM comparison_members WHERE comparison_id = ?", [cmp_id])
        DBInterface.execute(db, "DELETE FROM comparisons WHERE id = ?", [cmp_id])
        # Re-mint the comparison row at the original id so the dispatcher's
        # update branch lands on it.
        _premint_comparison!(db, cmp_id)

        HimalayaUI.rebuild_views_from_log!(db, cmp_id; entity_type="comparison")

        rebuilt_cmp = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT title, description, content_hash FROM comparisons WHERE id = ?", [cmp_id])))
        rebuilt_members = Tables.rowtable(DBInterface.execute(db,
            """SELECT id, exposure_id, display_order, snapshot
               FROM comparison_members WHERE comparison_id = ? ORDER BY id""", [cmp_id]))

        @test String(live_cmp.title) == String(rebuilt_cmp.title)
        @test String(live_cmp.description) == String(rebuilt_cmp.description)
        @test String(live_cmp.content_hash) == String(rebuilt_cmp.content_hash)
        @test length(live_members) == length(rebuilt_members)
        for (a, b) in zip(live_members, rebuilt_members)
            @test Int(a.exposure_id) == Int(b.exposure_id)
            @test Int(a.display_order) == Int(b.display_order)
            @test String(a.snapshot) == String(b.snapshot)
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# (Continuation of pre-existing tests below.)
# ─────────────────────────────────────────────────────────────────────────────

@testset "I2: apply_event!(::InTransaction, ...) rolls back with caller's tx on throw" begin
    mktempdir() do tmp
        db = HimalayaUI.open_db(joinpath(tmp, "test.db"))
        DBInterface.execute(db,
            "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e', '/p', '/d', '/a')")
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, name) VALUES (1, 'A1')")
        res = DBInterface.execute(db,
            "INSERT INTO exposures (sample_id, filename) VALUES (1, 'f')")
        exp_id = Int(DBInterface.lastrowid(res))

        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "op-roll"], UInt8[])
        try
            SQLite.transaction(db) do
                HimalayaUI.apply_event!(HimalayaUI.InTransaction(), db, req;
                    kind="peak_added", entity_type="exposure", entity_id=exp_id,
                    payload=Dict(:q => 3.0))
                error("intentional failure after apply_event!")
            end
        catch
        end
        # No row persisted — outer tx rolled back.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM user_actions WHERE client_op_id = 'op-roll'"))
        @test isempty(rows)
        crows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM peak_curations WHERE exposure_id = ?", [exp_id]))
        @test isempty(crows)
        SQLite.close(db)
    end
end

@testset "assignment_add round-trips through the log" begin
    mktempdir() do dir
        db  = HimalayaUI.open_db(joinpath(dir, "h.db"))
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        exp_id = HimalayaUI.create_experiment!(db; path="/x", data_dir="/x", analysis_dir="/x")
        s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id)
        e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id)
        DBInterface.execute(db, "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (10, ?, 'Pn3m', 0.1)", [e_id])

        HimalayaUI.apply_event!(db, req; kind="assignment_add",
            entity_type="exposure", entity_id=e_id, payload=Dict(:index_id => 10))

        # Live state: member present, state forced to 'indexed'.
        @test Set(Int(m.index_id) for m in Tables.rowtable(DBInterface.execute(db,
            "SELECT index_id FROM assignment_members WHERE exposure_id = ?", [e_id]))) == Set([10])
        @test String(Tables.rowtable(DBInterface.execute(db,
            "SELECT state FROM assignments WHERE exposure_id = ?", [e_id]))[1].state) == "indexed"

        # Wipe + rebuild from the log reproduces the member.
        DBInterface.execute(db, "DELETE FROM assignment_members WHERE exposure_id = ?", [e_id])
        HimalayaUI.rebuild_views_from_log!(db, e_id)
        @test Set(Int(m.index_id) for m in Tables.rowtable(DBInterface.execute(db,
            "SELECT index_id FROM assignment_members WHERE exposure_id = ?", [e_id]))) == Set([10])
    end
end
