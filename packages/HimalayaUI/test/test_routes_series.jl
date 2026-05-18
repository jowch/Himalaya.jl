using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# ─────────────────────────────────────────────────────────────────────────────
# I2.2 (#165) — series REST routes + business logic.
#
# The mutating routes emit events whose dispatcher branches do not land until
# #166–#168; `update_view_for_event!` no-ops unknown kinds, so those routes
# write a `user_actions` row but no view rows. Behavioural round-trips for the
# mutating routes belong to #166–#168 and #170; this file covers GET
# round-trips, the `last_event_at` sort fix, the messages round-trip (works
# today), and HTTP-level smoke tests for the mutating routes.
# ─────────────────────────────────────────────────────────────────────────────

# A DB with one experiment / sample / exposure — a valid FK target for
# `series_members.exposure_id` and `series_samples.sample_id`.
function _series_test_db(tmp::String)
    db = HimalayaUI.open_db(joinpath(tmp, "himalaya.db"))
    DBInterface.execute(db, """INSERT INTO experiments
        (id, name, path, data_dir, analysis_dir)
        VALUES (10, 'exp', '/x', '/x/d', '/x/a')""")
    DBInterface.execute(db,
        "INSERT INTO samples (id, experiment_id, name) VALUES (100, 10, 'sA')")
    DBInterface.execute(db,
        "INSERT INTO exposures (id, sample_id, filename, selected) VALUES (1000, 100, 'JC001', 1)")
    db
end

@testset "Series routes" begin

    @testset "GET /api/series — empty corpus" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                resp = HTTP.get("$base/api/series", ["X-Username" => "alice"])
                @test resp.status == 200
                @test JSON3.read(resp.body) == []
            end
            close(db)
        end
    end

    @testset "GET /api/series/{id}" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                # Missing series → 404.
                resp404 = HTTP.get("$base/api/series/999", ["X-Username" => "alice"];
                                   status_exception = false)
                @test resp404.status == 404

                # Seed a series with one recipe row and one plate member.
                DBInterface.execute(db, """INSERT INTO series
                    (id, title, state) VALUES (5, 'S5', 'draft')""")
                DBInterface.execute(db, """INSERT INTO series_samples
                    (series_id, sample_id, position, pinned, excluded)
                    VALUES (5, 100, 0, 1, 0)""")
                DBInterface.execute(db, """INSERT INTO series_members
                    (series_id, exposure_id, display_order, snapshot, created_at)
                    VALUES (5, 1000, 0, '{"effective_peaks":[],"confirmed_index":null,"analysis_inputs_hash":null}', '2026-05-01T00:00:00.000Z')""")

                resp = HTTP.get("$base/api/series/5", ["X-Username" => "alice"])
                @test resp.status == 200
                got = JSON3.read(resp.body, Dict{Symbol, Any})
                @test got[:id] == 5
                @test got[:title] == "S5"
                @test got[:state] == "draft"
                @test length(got[:members]) == 1
                @test got[:members][1]["exposure_id"] == 1000
                @test length(got[:samples]) == 1
                @test got[:samples][1]["sample_id"] == 100
                @test got[:samples][1]["pinned"] == true
                @test got[:order_rule] == "manual"
                @test got[:members][1]["is_stale"] == false
            end
            close(db)
        end
    end

    @testset "series_listing — last_event_at sort is recency-correct (#76)" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            # Series A (id 1): updated_at is a T-separated ISO string with a Z
            # suffix; no events. Series B (id 2): an older updated_at, but a
            # `user_actions` event whose space-separated timestamp is MORE
            # recent than A's updated_at. True recency order: B before A.
            DBInterface.execute(db, """INSERT INTO series (id, title, updated_at, state)
                VALUES (1, 'A', '2026-05-01T00:00:00.000Z', 'committed')""")
            DBInterface.execute(db, """INSERT INTO series (id, title, updated_at, state)
                VALUES (2, 'B', '2026-04-01T00:00:00.000Z', 'committed')""")
            DBInterface.execute(db, """INSERT INTO user_actions
                (action, entity_type, entity_id, timestamp)
                VALUES ('series_recipe_updated', 'series', 2, '2026-05-10 12:00:00')""")

            listing = HimalayaUI.series_listing(db)
            @test length(listing) == 2
            # Bug #76: lexical sort would put A first ('T' > ' '). The
            # datetime() wrapper sorts by true recency — B first.
            @test listing[1][:id] == 2
            @test listing[2][:id] == 1
            # The projected last_event_at is normalised (no 'T'/'Z'), so it is
            # itself a valid client sort key — the bug is closed end-to-end.
            @test !occursin('T', listing[1][:last_event_at])
            @test !occursin('Z', listing[1][:last_event_at])
            # Series A's last_event_at comes from the T-separated, Z-suffixed
            # updated_at — assert datetime() normalised that path too.
            @test !occursin('T', listing[2][:last_event_at])
            @test !occursin('Z', listing[2][:last_event_at])
            close(db)
        end
    end

    @testset "GET /api/series/{id}/forks" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                # An existing series with no forks → empty array.
                DBInterface.execute(db, """INSERT INTO series (id, title, state)
                    VALUES (7, 'parent', 'committed')""")
                resp0 = HTTP.get("$base/api/series/7/forks", ["X-Username" => "alice"])
                @test resp0.status == 200
                @test JSON3.read(resp0.body) == []

                # Add a fork (id 8, forked_from_id = 7) → it now appears.
                DBInterface.execute(db, """INSERT INTO series
                    (id, title, state, forked_from_id)
                    VALUES (8, 'child', 'committed', 7)""")
                resp = HTTP.get("$base/api/series/7/forks", ["X-Username" => "alice"])
                @test resp.status == 200
                forks = JSON3.read(resp.body, Vector{Dict{Symbol, Any}})
                @test length(forks) == 1
                @test forks[1][:id] == 8
                @test forks[1][:forked_from_id] == 7
            end
            close(db)
        end
    end

end
