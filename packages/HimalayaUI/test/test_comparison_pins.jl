using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# ─────────────────────────────────────────────────────────────────────────────
# Phase 13 Task 13.2 — comparison pins.
#
# Per-user pinned comparisons surface at the top of the sidebar. The backend
# is a small `comparison_pins(user_id, comparison_id, pinned_at)` table with
# composite PK; pin/unpin are trivial idempotent state toggles, so the routes
# do NOT wrap in `with_idempotency`. The routes return a small JSON payload
# echoing `{comparison_id, pinned}`.
#
# Helpers (`_premint_cmp!`, `_setup_analyzed_exposure`) are sourced from
# test_comparisons.jl which loaded earlier in runtests.jl.
# ─────────────────────────────────────────────────────────────────────────────

@testset "comparison pins" begin

    @testset "schema: comparison_pins exists with composite PK + index" begin
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
            try
                # The migration runs on open_db; the table should exist.
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT 1 FROM sqlite_master WHERE type='table' AND name='comparison_pins'"))
                @test !isempty(rows)
                # Index for per-user listing.
                idx = Tables.rowtable(DBInterface.execute(db,
                    "SELECT 1 FROM sqlite_master WHERE type='index' AND name='idx_comparison_pins_by_user'"))
                @test !isempty(idx)
            finally
                close(db)
            end
        end
    end

    @testset "POST .../pin: 401 when X-Username missing" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 1)
            with_test_server(ctx.db) do port, base
                r = HTTP.post("$base/api/comparisons/1/pin"; status_exception=false)
                @test r.status == 401
            end
        end
    end

    @testset "POST .../pin: 404 when comparison missing" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.post("$base/api/comparisons/999/pin",
                              ["X-Username" => "alice"]; status_exception=false)
                @test r.status == 404
            end
        end
    end

    @testset "POST .../pin → DELETE .../pin round-trip persists in table" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 1)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET content_hash = 'h1' WHERE id = 1")
            with_test_server(ctx.db) do port, base
                r1 = HTTP.post("$base/api/comparisons/1/pin",
                              ["X-Username" => "alice"])
                @test r1.status == 200
                body1 = JSON3.read(String(r1.body))
                @test body1.pinned == true
                @test body1.comparison_id == 1

                # Verify in DB.
                rows = Tables.rowtable(DBInterface.execute(ctx.db,
                    "SELECT comparison_id FROM comparison_pins"))
                @test length(rows) == 1
                @test Int(rows[1].comparison_id) == 1

                # Unpin returns pinned=false and removes the row.
                r2 = HTTP.delete("$base/api/comparisons/1/pin",
                                 ["X-Username" => "alice"])
                @test r2.status == 200
                body2 = JSON3.read(String(r2.body))
                @test body2.pinned == false

                rows2 = Tables.rowtable(DBInterface.execute(ctx.db,
                    "SELECT comparison_id FROM comparison_pins"))
                @test isempty(rows2)
            end
        end
    end

    @testset "POST .../pin twice: idempotent (same row, refreshed timestamp)" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 1)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET content_hash = 'h1' WHERE id = 1")
            with_test_server(ctx.db) do port, base
                HTTP.post("$base/api/comparisons/1/pin", ["X-Username" => "alice"])
                # Second pin — must not duplicate the row.
                HTTP.post("$base/api/comparisons/1/pin", ["X-Username" => "alice"])
                rows = Tables.rowtable(DBInterface.execute(ctx.db,
                    "SELECT COUNT(*) AS n FROM comparison_pins"))
                @test Int(rows[1].n) == 1
            end
        end
    end

    @testset "DELETE .../pin: idempotent — never-pinned comparison still 200" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 1)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET content_hash = 'h1' WHERE id = 1")
            with_test_server(ctx.db) do port, base
                r = HTTP.delete("$base/api/comparisons/1/pin",
                                ["X-Username" => "alice"])
                @test r.status == 200
                body = JSON3.read(String(r.body))
                @test body.pinned == false
            end
        end
    end

    @testset "GET /api/users/me/comparison-pins: 401 without X-Username" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/users/me/comparison-pins";
                             status_exception=false)
                @test r.status == 401
            end
        end
    end

    @testset "GET /api/users/me/comparison-pins: empty array for new user" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/users/me/comparison-pins",
                             ["X-Username" => "alice"])
                @test r.status == 200
                @test JSON3.read(String(r.body)) == []
            end
        end
    end

    @testset "GET /api/users/me/comparison-pins: most-recently-pinned first" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 1)
            _premint_cmp!(ctx.db, 2)
            _premint_cmp!(ctx.db, 3)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET content_hash = 'h' WHERE id IN (1,2,3)")
            # Manually set distinct pinned_at timestamps so ordering is
            # deterministic regardless of clock resolution. The CURRENT_TIMESTAMP
            # default in the route uses second-resolution; rapid-fire calls can
            # share a timestamp and tie-break by comparison_id DESC instead.
            user_id = HimalayaUI.get_or_create_user!(ctx.db, "alice")
            DBInterface.execute(ctx.db, """
                INSERT INTO comparison_pins (user_id, comparison_id, pinned_at) VALUES
                (?, 2, '2026-05-01T00:00:00'),
                (?, 3, '2026-05-01T00:00:01'),
                (?, 1, '2026-05-01T00:00:02')
            """, [user_id, user_id, user_id])
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/users/me/comparison-pins",
                             ["X-Username" => "alice"])
                @test r.status == 200
                ids = JSON3.read(String(r.body))
                # Most-recent first: 1, 3, 2.
                @test ids == [1, 3, 2]
            end
        end
    end

    @testset "GET /api/users/me/comparison-pins: per-user isolation" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 1)
            _premint_cmp!(ctx.db, 2)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET content_hash = 'h' WHERE id IN (1,2)")
            with_test_server(ctx.db) do port, base
                HTTP.post("$base/api/comparisons/1/pin", ["X-Username" => "alice"])
                HTTP.post("$base/api/comparisons/2/pin", ["X-Username" => "bob"])
                r_a = HTTP.get("$base/api/users/me/comparison-pins",
                               ["X-Username" => "alice"])
                r_b = HTTP.get("$base/api/users/me/comparison-pins",
                               ["X-Username" => "bob"])
                @test JSON3.read(String(r_a.body)) == [1]
                @test JSON3.read(String(r_b.body)) == [2]
            end
        end
    end

    # ─── Event-log round-trip (Phase 4 follow-up) ─────────────────────────────
    #
    # Phase 13 originally landed pin/unpin as direct SQL writes — no
    # event-log row, no SSE broadcast. The follow-up routes everything
    # through `apply_event!` so cross-tab fanout works and the durable
    # history is queryable. These tests pin the round-trip property: live
    # writes via `apply_event!` produce the same view-table state as
    # `rebuild_views_from_log!` rebuilding from `user_actions` alone.

    @testset "comparison_pinned: rebuild_views_from_log! matches live state" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 1)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET content_hash = 'h1' WHERE id = 1")

            user_id = HimalayaUI.get_or_create_user!(ctx.db, "alice")
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            HimalayaUI.apply_event!(ctx.db, req;
                kind = "comparison_pinned",
                entity_type = "user",
                entity_id = user_id,
                payload = Dict(:comparison_id => 1))

            before = Tables.rowtable(DBInterface.execute(ctx.db, """
                SELECT user_id, comparison_id FROM comparison_pins
                ORDER BY user_id, comparison_id
            """))
            @test length(before) == 1
            @test Int(before[1].comparison_id) == 1

            # Wipe + replay from log.
            DBInterface.execute(ctx.db, "DELETE FROM comparison_pins")
            HimalayaUI.rebuild_views_from_log!(ctx.db, user_id; entity_type = "user")

            after = Tables.rowtable(DBInterface.execute(ctx.db, """
                SELECT user_id, comparison_id FROM comparison_pins
                ORDER BY user_id, comparison_id
            """))
            @test length(after) == length(before)
            @test Int(after[1].comparison_id) == Int(before[1].comparison_id)
            @test Int(after[1].user_id)       == Int(before[1].user_id)
        end
    end

    @testset "comparison_unpinned: rebuild_views_from_log! matches live state" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 1)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET content_hash = 'h1' WHERE id = 1")
            user_id = HimalayaUI.get_or_create_user!(ctx.db, "alice")
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])

            # Pin then unpin via apply_event!.
            HimalayaUI.apply_event!(ctx.db, req;
                kind = "comparison_pinned",
                entity_type = "user",
                entity_id = user_id,
                payload = Dict(:comparison_id => 1))
            HimalayaUI.apply_event!(ctx.db, req;
                kind = "comparison_unpinned",
                entity_type = "user",
                entity_id = user_id,
                payload = Dict(:comparison_id => 1))

            before = Tables.rowtable(DBInterface.execute(ctx.db,
                "SELECT 1 FROM comparison_pins"))
            @test isempty(before)

            # Reset views, rebuild from log — pin then unpin → no rows.
            DBInterface.execute(ctx.db, "DELETE FROM comparison_pins")
            HimalayaUI.rebuild_views_from_log!(ctx.db, user_id; entity_type = "user")

            after = Tables.rowtable(DBInterface.execute(ctx.db,
                "SELECT 1 FROM comparison_pins"))
            @test isempty(after)
        end
    end

    @testset "POST .../pin writes a comparison_pinned event row" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 1)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET content_hash = 'h1' WHERE id = 1")
            with_test_server(ctx.db) do port, base
                r = HTTP.post("$base/api/comparisons/1/pin",
                              ["X-Username" => "alice"])
                @test r.status == 200
                # One event row, kind=comparison_pinned, entity_type=user,
                # entity_id=user_id, payload contains comparison_id=1.
                rows = Tables.rowtable(DBInterface.execute(ctx.db, """
                    SELECT action, entity_type, entity_id, payload FROM user_actions
                    WHERE action = 'comparison_pinned'
                """))
                @test length(rows) == 1
                @test String(rows[1].entity_type) == "user"
                payload = JSON3.read(String(rows[1].payload))
                @test Int(payload.comparison_id) == 1
            end
        end
    end

    @testset "DELETE .../pin writes a comparison_unpinned event row" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 1)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET content_hash = 'h1' WHERE id = 1")
            with_test_server(ctx.db) do port, base
                HTTP.post("$base/api/comparisons/1/pin", ["X-Username" => "alice"])
                r = HTTP.delete("$base/api/comparisons/1/pin",
                                ["X-Username" => "alice"])
                @test r.status == 200
                rows = Tables.rowtable(DBInterface.execute(ctx.db, """
                    SELECT entity_type, payload FROM user_actions
                    WHERE action = 'comparison_unpinned'
                """))
                @test length(rows) == 1
                @test String(rows[1].entity_type) == "user"
                payload = JSON3.read(String(rows[1].payload))
                @test Int(payload.comparison_id) == 1
            end
        end
    end

    @testset "FK ON DELETE CASCADE: pin disappears when comparison deleted" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 1)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET content_hash = 'h1' WHERE id = 1")
            with_test_server(ctx.db) do port, base
                HTTP.post("$base/api/comparisons/1/pin", ["X-Username" => "alice"])
                # Hard-delete the comparison directly (bypassing route auth) to
                # exercise the cascade.
                DBInterface.execute(ctx.db, "DELETE FROM comparisons WHERE id = 1")
                rows = Tables.rowtable(DBInterface.execute(ctx.db,
                    "SELECT 1 FROM comparison_pins"))
                @test isempty(rows)
            end
        end
    end

end
