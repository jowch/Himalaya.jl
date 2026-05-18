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

    @testset "content-hash + existence helpers" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            # series_exists / current_series_content_hash on a missing id.
            @test HimalayaUI.series_exists(db, 404) == false
            @test HimalayaUI.current_series_content_hash(db, 404) === nothing

            # A draft series: it exists, but its content_hash is NULL — the
            # reason series_exists must be a separate probe (spec §4).
            DBInterface.execute(db, """INSERT INTO series (id, title, state)
                VALUES (3, 'draft-s', 'draft')""")
            @test HimalayaUI.series_exists(db, 3) == true
            @test HimalayaUI.current_series_content_hash(db, 3) === nothing

            # compute_series_content_hash is plate-only: adding a recipe row
            # (series_samples) must NOT change the hash; adding a plate member
            # (series_members) must.
            h_empty = HimalayaUI.compute_series_content_hash(db, 3)
            DBInterface.execute(db, """INSERT INTO series_samples
                (series_id, sample_id, position) VALUES (3, 100, 0)""")
            @test HimalayaUI.compute_series_content_hash(db, 3) == h_empty
            DBInterface.execute(db, """INSERT INTO series_members
                (series_id, exposure_id, display_order, snapshot, created_at)
                VALUES (3, 1000, 0, '{"effective_peaks":[],"confirmed_index":null,"analysis_inputs_hash":null}', '2026-05-01T00:00:00.000Z')""")
            @test HimalayaUI.compute_series_content_hash(db, 3) != h_empty

            close(db)
        end
    end

    @testset "POST /api/series — create draft (smoke)" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                # Missing title → 400 (uncached validation error).
                resp400 = HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:samples => []));
                    status_exception = false)
                @test resp400.status == 400

                # samples present but not an array → 400.
                resp400b = HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:title => "t", :samples => "nope"));
                    status_exception = false)
                @test resp400b.status == 400

                # Well-formed create → 201; a series row + a user_actions row.
                resp = HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(
                        :title => "My series",
                        :order_rule => "ascending",
                        :samples => [Dict(:sample_id => 100, :position => 0,
                                          :pinned => false, :excluded => false)])))
                @test resp.status == 201
                created = JSON3.read(resp.body, Dict{Symbol, Any})
                new_id = created[:id]
                @test new_id isa Integer
                # Degenerate until #166: the dispatcher no-ops, so the body is
                # the placeholder projection — empty members, empty samples.
                @test created[:members] == []
                @test created[:samples] == []
                # The durable event row IS written.
                ev = Tables.rowtable(DBInterface.execute(db,
                    "SELECT action FROM user_actions WHERE entity_type='series' AND entity_id=?",
                    [new_id]))
                @test length(ev) == 1
                @test ev[1].action == "series_created"

                # An empty samples array is a valid draft — the deliberate
                # departure from comparisons, which rejects empty members.
                respEmpty = HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:title => "empty draft", :samples => [])))
                @test respEmpty.status == 201

                # A samples entry missing sample_id → uncached 400.
                resp400c = HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:title => "t", :samples => [Dict(:position => 0)]));
                    status_exception = false)
                @test resp400c.status == 400
            end
            close(db)
        end
    end

    @testset "PATCH /api/series/{id} — recipe edit (smoke)" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                # Missing series → 404.
                resp404 = HTTP.patch("$base/api/series/999",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:samples => []));
                    status_exception = false)
                @test resp404.status == 404

                # samples not an array → 400.
                DBInterface.execute(db, """INSERT INTO series (id, title, state)
                    VALUES (12, 'edit-me', 'draft')""")
                resp400 = HTTP.patch("$base/api/series/12",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:samples => "nope"));
                    status_exception = false)
                @test resp400.status == 400

                # A samples entry missing sample_id → uncached 400.
                resp400b = HTTP.patch("$base/api/series/12",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:samples => [Dict(:position => 0)]));
                    status_exception = false)
                @test resp400b.status == 400

                # Well-formed recipe edit → 200; a user_actions row written.
                resp = HTTP.patch("$base/api/series/12",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(
                        :order_rule => "descending",
                        :samples => [Dict(:sample_id => 100, :position => 0)])))
                @test resp.status == 200
                # Degenerate until #166: the dispatcher no-ops, so the recipe
                # is not applied — only assert the body is series 12's projection.
                @test JSON3.read(resp.body, Dict{Symbol, Any})[:id] == 12
                ev = Tables.rowtable(DBInterface.execute(db,
                    "SELECT action FROM user_actions WHERE entity_type='series' AND entity_id=12"))
                @test length(ev) == 1
                @test ev[1].action == "series_recipe_updated"
            end
            close(db)
        end
    end

    @testset "POST /api/series/{id}/commit (smoke + no gate + 409)" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                _member = Dict(:exposure_id => 1000, :display_order => 0,
                               :snapshot => Dict(:effective_peaks => Any[],
                                                 :confirmed_index => nothing,
                                                 :analysis_inputs_hash => nothing))

                # Missing series → 404.
                resp404 = HTTP.post("$base/api/series/999/commit",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:members => [_member]));
                    status_exception = false)
                @test resp404.status == 404

                # Seed a committed series authored by alice (id 1), with a
                # stored content_hash so the 409 path can be exercised.
                DBInterface.execute(db, "INSERT INTO users (id, username) VALUES (1, 'alice')")
                DBInterface.execute(db, """INSERT INTO series
                    (id, title, state, created_by, content_hash)
                    VALUES (20, 'committed-s', 'committed', 1, 'sha256:deadbeef')""")

                # No is_author gate: bob (not the author) commits → NOT 403.
                resp = HTTP.post("$base/api/series/20/commit",
                    ["X-Username" => "bob", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:members => [_member]));
                    status_exception = false)
                @test resp.status != 403
                @test resp.status == 200
                ev = Tables.rowtable(DBInterface.execute(db,
                    "SELECT action FROM user_actions WHERE entity_type='series' AND entity_id=20"))
                @test any(r -> r.action == "series_plate_committed", ev)

                # Conflict: a wrong expected_content_hash → 409.
                resp409 = HTTP.post("$base/api/series/20/commit",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:members => [_member],
                                     :expected_content_hash => "sha256:WRONG"));
                    status_exception = false)
                @test resp409.status == 409
                conflict = JSON3.read(resp409.body, Dict{Symbol, Any})
                @test conflict[:error] == "conflict"
                @test conflict[:current_hash] == "sha256:deadbeef"
            end
            close(db)
        end
    end

end
