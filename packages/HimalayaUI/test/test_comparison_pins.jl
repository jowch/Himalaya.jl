using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# ─────────────────────────────────────────────────────────────────────────────
# Phase 13 Task 13.2 — comparison pins.
#
# I3.6 (#177): the Compare page + its `/api/comparisons/:id/pin`,
# `/unpin`, and `/api/users/me/comparison-pins` routes are retired with
# `routes_comparisons.jl`. The route-level tests (401/404/round-trip/idempotent/
# listing/event-row/per-user-isolation) went with them. What stays here is the
# KEPT-FOREVER event-replay machinery: the `comparison_pins` table, the
# `comparison_pinned` / `comparison_unpinned` dispatcher branches, and the
# `rebuild_views_from_log!` fold — driven directly via `apply_event!` (no HTTP),
# so a migrated/historical pin event still replays after the routes are gone.
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

    @testset "FK ON DELETE CASCADE: pin disappears when comparison deleted" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            _premint_cmp!(ctx.db, 1)
            DBInterface.execute(ctx.db,
                "UPDATE comparisons SET content_hash = 'h1' WHERE id = 1")
            user_id = HimalayaUI.get_or_create_user!(ctx.db, "alice")
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            # Pin via apply_event! (the route is retired; the dispatcher fold is
            # what writes the comparison_pins row).
            HimalayaUI.apply_event!(ctx.db, req;
                kind = "comparison_pinned",
                entity_type = "user",
                entity_id = user_id,
                payload = Dict(:comparison_id => 1))
            # Hard-delete the comparison directly to exercise the schema's
            # ON DELETE CASCADE on comparison_pins.
            DBInterface.execute(ctx.db, "DELETE FROM comparisons WHERE id = 1")
            rows = Tables.rowtable(DBInterface.execute(ctx.db,
                "SELECT 1 FROM comparison_pins"))
            @test isempty(rows)
        end
    end

end
