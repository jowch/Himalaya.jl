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

# Capture the SSE frames of one kind broadcast while `f` runs — the
# in-process subscriber pattern (test/AGENTS.md). `do`-block form: `f` is the
# FIRST argument, so `_capture_series_sse("k") do … end` ⇒ (f, "k").
function _capture_series_sse(f::Function, kind::String)
    pending = Channel{String}(64)
    sub = (pending = pending,)
    lock(HimalayaUI.SSE_LOCK) do
        push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
    end
    try
        f()
        sleep(0.3)   # let the post-commit broadcast queue flush
    finally
        lock(HimalayaUI.SSE_LOCK) do
            filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
        end
        close(pending)
    end
    [fr for fr in pending
        if !startswith(fr, ":") && occursin("\"kind\":\"$kind\"", fr)]
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
                # #166: dispatcher now writes state='draft' + recipe rows.
                @test created[:state] == "draft"
                @test created[:members] == []
                @test length(created[:samples]) == 1
                # The durable event row IS written.
                ev = Tables.rowtable(DBInterface.execute(db,
                    "SELECT action FROM user_actions WHERE entity_type='series' AND entity_id=?",
                    [new_id]))
                @test length(ev) == 1
                @test ev[1].action == "series_created"

                # Lock the _series_sample_payload shape in the durable event
                # payload — #166 builds its dispatcher against this shape.
                evrow = Tables.rowtable(DBInterface.execute(db,
                    "SELECT payload FROM user_actions WHERE entity_type='series' AND entity_id=?",
                    [new_id]))[1]
                evpayload = JSON3.read(String(evrow.payload))
                @test length(evpayload.samples) == 1
                @test evpayload.samples[1].sample_id == 100
                @test evpayload.samples[1].position == 0
                @test evpayload.samples[1].pinned == false
                @test evpayload.samples[1].excluded == false

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
                # Capture the hash the dispatcher just wrote (plate-based, not the seed).
                committed_hash = JSON3.read(resp.body, Dict{Symbol, Any})[:content_hash]

                # Conflict: a wrong expected_content_hash → 409.
                resp409 = HTTP.post("$base/api/series/20/commit",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:members => [_member],
                                     :expected_content_hash => "sha256:WRONG"));
                    status_exception = false)
                @test resp409.status == 409
                conflict = JSON3.read(resp409.body, Dict{Symbol, Any})
                @test conflict[:error] == "conflict"
                @test conflict[:current_hash] == committed_hash
            end
            close(db)
        end
    end

    @testset "DELETE /api/series/{id} (smoke + no gate)" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                # Missing series → 404.
                resp404 = HTTP.delete("$base/api/series/999", ["X-Username" => "alice"];
                                      status_exception = false)
                @test resp404.status == 404

                # Series authored by alice; bob (not the author) deletes it —
                # no is_author gate, so NOT 403.
                DBInterface.execute(db, "INSERT INTO users (id, username) VALUES (1, 'alice')")
                DBInterface.execute(db, """INSERT INTO series
                    (id, title, state, created_by) VALUES (30, 'doomed', 'committed', 1)""")
                resp = HTTP.delete("$base/api/series/30", ["X-Username" => "bob"];
                                   status_exception = false)
                @test resp.status != 403
                @test resp.status == 200
                deleted = JSON3.read(resp.body, Dict{Symbol, Any})
                @test deleted[:id] == 30
                @test deleted[:deleted] == true
                @test deleted[:event_id] isa Integer && deleted[:event_id] > 0
                ev = Tables.rowtable(DBInterface.execute(db,
                    "SELECT action FROM user_actions WHERE entity_type='series' AND entity_id=30"))
                @test any(r -> r.action == "series_deleted", ev)
            end
            close(db)
        end
    end

    @testset "series messages — full round-trip" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                DBInterface.execute(db, """INSERT INTO series (id, title, state)
                    VALUES (40, 'chatty', 'committed')""")

                # Empty thread.
                resp0 = HTTP.get("$base/api/series/40/messages", ["X-Username" => "alice"])
                @test resp0.status == 200
                @test JSON3.read(resp0.body) == []

                # POST without X-Username → 401.
                resp401 = HTTP.post("$base/api/series/40/messages",
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:body => "hi"));
                    status_exception = false)
                @test resp401.status == 401

                # POST with an empty body → 400.
                resp400 = HTTP.post("$base/api/series/40/messages",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:body => "   "));
                    status_exception = false)
                @test resp400.status == 400

                # POST a real message → 201, then GET sees it.
                resp = HTTP.post("$base/api/series/40/messages",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:body => "first post")))
                @test resp.status == 201
                msg = JSON3.read(resp.body, Dict{Symbol, Any})
                @test msg[:body] == "first post"
                @test msg[:series_id] == 40

                resp2 = HTTP.get("$base/api/series/40/messages", ["X-Username" => "alice"])
                thread = JSON3.read(resp2.body, Vector{Dict{Symbol, Any}})
                @test length(thread) == 1
                @test thread[1][:body] == "first post"
                @test thread[1][:author] == "alice"
            end
            close(db)
        end
    end

    @testset "series pins (smoke)" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                # series-pins without X-Username → 401.
                resp401 = HTTP.get("$base/api/users/me/series-pins";
                                   status_exception = false)
                @test resp401.status == 401

                # Fresh user → empty pin list.
                resp0 = HTTP.get("$base/api/users/me/series-pins", ["X-Username" => "alice"])
                @test resp0.status == 200
                @test JSON3.read(resp0.body) == []

                # Pin a missing series → 404.
                resp404 = HTTP.post("$base/api/series/999/pin", ["X-Username" => "alice"];
                                    status_exception = false)
                @test resp404.status == 404

                # Pin an existing series → 200, a user_actions row written.
                DBInterface.execute(db, """INSERT INTO series (id, title, state)
                    VALUES (50, 'pin-me', 'committed')""")
                respPin = HTTP.post("$base/api/series/50/pin", ["X-Username" => "alice"])
                @test respPin.status == 200
                pinned = JSON3.read(respPin.body, Dict{Symbol, Any})
                @test pinned[:series_id] == 50
                @test pinned[:pinned] == true
                evP = Tables.rowtable(DBInterface.execute(db,
                    "SELECT action FROM user_actions WHERE entity_type='user'"))
                @test any(r -> r.action == "series_pinned", evP)

                # Unpin → 200.
                respUnpin = HTTP.delete("$base/api/series/50/pin", ["X-Username" => "alice"])
                @test respUnpin.status == 200
                unpinned = JSON3.read(respUnpin.body, Dict{Symbol, Any})
                @test unpinned[:pinned] == false
                evU = Tables.rowtable(DBInterface.execute(db,
                    "SELECT action FROM user_actions WHERE entity_type='user'"))
                @test any(r -> r.action == "series_unpinned", evU)
            end
            close(db)
        end
    end

    @testset "POST /api/series — series_created writes the recipe + draft state" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                body = JSON3.write(Dict(
                    :title   => "Heat ramp",
                    :samples => [Dict(:sample_id => 100, :position => 0, :pinned => true)],
                ))
                resp = HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"], body)
                @test resp.status == 201
                got = JSON3.read(resp.body, Dict{Symbol, Any})
                sid = got[:id]
                # The dispatcher upserts state='draft' and leaves content_hash NULL.
                @test got[:state] == "draft"
                @test got[:content_hash] == ""           # fetch_series_with_plate maps NULL → ""
                @test got[:title] == "Heat ramp"
                # The recipe snapshot landed in series_samples.
                @test length(got[:samples]) == 1
                @test got[:samples][1]["sample_id"] == 100
                @test got[:samples][1]["pinned"] == true
                @test isempty(got[:members])             # series_created carries zero members

                # rebuild_views_from_log! round-trip: empty the view rows, re-fold.
                DBInterface.execute(db, "DELETE FROM series_samples WHERE series_id = ?", [sid])
                DBInterface.execute(db, "DELETE FROM series WHERE id = ?", [sid])
                HimalayaUI.rebuild_views_from_log!(db, sid; entity_type = "series")
                refold = HimalayaUI.fetch_series_with_plate(db, sid)
                @test refold !== nothing
                @test refold[:state] == "draft"
                @test length(refold[:samples]) == 1
                @test refold[:samples][1][:sample_id] == 100

                # SSE layer: a second create, observed through the in-process
                # subscriber, must broadcast exactly one series_created frame
                # carrying entity_type='series'.
                frames = _capture_series_sse("series_created") do
                    HTTP.post("$base/api/series",
                        ["X-Username" => "alice", "Content-Type" => "application/json"],
                        JSON3.write(Dict(:title => "Frame check",
                            :samples => [Dict(:sample_id => 100, :position => 0)])))
                end
                @test length(frames) == 1
                @test occursin("\"entity_type\":\"series\"", frames[1])
            end
            close(db)
        end
    end

    @testset "PATCH /api/series/{id} — series_recipe_updated pure-replaces the recipe" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            DBInterface.execute(db,
                "INSERT INTO samples (id, experiment_id, name) VALUES (101, 10, 'sB')")
            with_test_server(db) do port, base
                # Create a draft with one recipe row.
                createBody = JSON3.write(Dict(
                    :title   => "Ramp",
                    :samples => [Dict(:sample_id => 100, :position => 0)],
                ))
                created = JSON3.read(HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    createBody).body, Dict{Symbol, Any})
                sid = created[:id]

                # PATCH with a completely different recipe — pure-replace, not a diff.
                patchBody = JSON3.write(Dict(
                    :ordering_variable => "temperature",
                    :order_rule        => "ascending",
                    :samples           => [
                        Dict(:sample_id => 101, :position => 0),
                        Dict(:sample_id => 100, :position => 1),
                    ],
                ))
                resp = HTTP.patch("$base/api/series/$sid",
                    ["X-Username" => "alice", "Content-Type" => "application/json"], patchBody)
                @test resp.status == 200
                got = JSON3.read(resp.body, Dict{Symbol, Any})
                @test got[:ordering_variable] == "temperature"
                @test got[:order_rule] == "ascending"
                @test length(got[:samples]) == 2
                @test got[:samples][1]["sample_id"] == 101   # ordered by position
                @test got[:samples][2]["sample_id"] == 100
                @test got[:state] == "draft"                  # recipe edit never commits
                @test got[:content_hash] == ""                # recipe edit never hashes

                # rebuild_views_from_log! round-trip: fold series_created THEN
                # series_recipe_updated from empty; final state == post-PATCH.
                DBInterface.execute(db, "DELETE FROM series_samples WHERE series_id = ?", [sid])
                DBInterface.execute(db, "DELETE FROM series WHERE id = ?", [sid])
                HimalayaUI.rebuild_views_from_log!(db, sid; entity_type = "series")
                refold = HimalayaUI.fetch_series_with_plate(db, sid)
                @test length(refold[:samples]) == 2
                @test refold[:samples][1][:sample_id] == 101
                @test refold[:ordering_variable] == "temperature"
            end
            close(db)
        end
    end

    @testset "DELETE /api/series/{id} — series_deleted cascades + folds to absent" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                createBody = JSON3.write(Dict(
                    :title   => "Doomed",
                    :samples => [Dict(:sample_id => 100, :position => 0)],
                ))
                sid = JSON3.read(HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    createBody).body, Dict{Symbol, Any})[:id]

                resp = HTTP.delete("$base/api/series/$sid", ["X-Username" => "alice"])
                @test resp.status == 200
                # Four-table cascade: the series row and its recipe rows are gone.
                @test HimalayaUI.fetch_series_with_plate(db, sid) === nothing
                remaining = Tables.rowtable(DBInterface.execute(db,
                    "SELECT COUNT(*) AS n FROM series_samples WHERE series_id = ?", [sid]))
                @test remaining[1].n == 0

                # rebuild_views_from_log! round-trip: fold series_created THEN
                # series_deleted from empty → series stays absent.
                HimalayaUI.rebuild_views_from_log!(db, sid; entity_type = "series")
                @test HimalayaUI.fetch_series_with_plate(db, sid) === nothing
            end
            close(db)
        end
    end

    @testset "POST /api/series/{id}/commit — series_plate_committed commits the plate" begin
        mktempdir() do tmp
            db = _series_test_db(tmp)
            with_test_server(db) do port, base
                sid = JSON3.read(HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:title => "Ramp",
                                     :samples => [Dict(:sample_id => 100, :position => 0)]))
                    ).body, Dict{Symbol, Any})[:id]

                # Pass an explicit snapshot so the route's _series_member_payload
                # uses it directly (no dependence on compute_member_snapshot
                # succeeding against the unanalyzed exposure fixture).
                commitBody = JSON3.write(Dict(:members => [
                    Dict(:exposure_id => 1000, :display_order => 0,
                         :snapshot => Dict(:effective_peaks => [],
                                           :confirmed_index => nothing,
                                           :analysis_inputs_hash => nothing)),
                ]))
                resp = HTTP.post("$base/api/series/$sid/commit",
                    ["X-Username" => "alice", "Content-Type" => "application/json"], commitBody)
                @test resp.status == 200
                got = JSON3.read(resp.body, Dict{Symbol, Any})
                @test got[:state] == "committed"
                @test startswith(got[:content_hash], "sha256:")   # hashed from the plate
                @test length(got[:members]) == 1
                @test got[:members][1]["exposure_id"] == 1000
                @test length(got[:samples]) == 1                  # recipe untouched

                # rebuild_views_from_log! round-trip from empty.
                DBInterface.execute(db, "DELETE FROM series_members WHERE series_id = ?", [sid])
                DBInterface.execute(db, "DELETE FROM series_samples WHERE series_id = ?", [sid])
                DBInterface.execute(db, "DELETE FROM series WHERE id = ?", [sid])
                HimalayaUI.rebuild_views_from_log!(db, sid; entity_type = "series")
                refold = HimalayaUI.fetch_series_with_plate(db, sid)
                @test refold[:state] == "committed"
                @test length(refold[:members]) == 1
                @test refold[:content_hash] == got[:content_hash]   # hash is fold-stable

                # SSE layer: series_plate_committed is the one series event
                # carrying a post_state envelope. Create + commit a fresh
                # series through the in-process subscriber and assert the
                # frame carries entity_type='series' and a post_state field.
                sid2 = JSON3.read(HTTP.post("$base/api/series",
                    ["X-Username" => "alice", "Content-Type" => "application/json"],
                    JSON3.write(Dict(:title => "Frame check",
                        :samples => [Dict(:sample_id => 100, :position => 0)]))).body,
                    Dict{Symbol, Any})[:id]
                frames = _capture_series_sse("series_plate_committed") do
                    HTTP.post("$base/api/series/$sid2/commit",
                        ["X-Username" => "alice", "Content-Type" => "application/json"],
                        JSON3.write(Dict(:members => [
                            Dict(:exposure_id => 1000, :display_order => 0,
                                 :snapshot => Dict(:effective_peaks => [],
                                                   :confirmed_index => nothing,
                                                   :analysis_inputs_hash => nothing))])))
                end
                @test length(frames) == 1
                @test occursin("\"entity_type\":\"series\"", frames[1])
                @test occursin("\"post_state\"", frames[1])
            end
            close(db)
        end
    end

end
