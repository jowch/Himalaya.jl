using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# Idempotency replay invariant (Tier 1.4).
#
# Per-route, posting the same X-Client-Op-Id twice MUST produce:
#   1. Byte-identical HTTP response bodies (already checked at unit level
#      in test_idempotency.jl, but pinned end-to-end through real routes).
#   2. Exactly ONE durable row in user_actions (the second call is a cache
#      replay — no new event row).
#   3. Exactly ONE SSE frame fanout (the second call must NOT broadcast a
#      duplicate frame; subscribers depend on this for own-op deferred
#      resolution to be unique-keyed).
#
# Existing test_idempotency.jl unit tests pin (1) at the with_idempotency
# helper level. Existing per-route tests cover (2) implicitly. Pinning all
# three together at the route level — across the curation, group, and
# message kinds — closes the contract.
#
# What's NOT tested here: concurrent replay (test_idempotency.jl already
# covers OP_LOCKS serialization).

# Helper — count user_actions rows + filter by kind.
function _count_actions(db, kind)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM user_actions WHERE action = ?", [kind]))
    Int(rows[1].c)
end

# Helper — register a direct in-process SSE subscriber (no HTTP), run the
# closure, then drain captured frames. Subscribers in HimalayaUI.SSE_SUBSCRIBERS
# receive a JSON-encoded `event: curation\ndata: {...}\n\n` string per
# `broadcast_event!` call. Filtering by `kind_filter` (the event kind from
# the JSON) lets us count only the frames the test cares about.
#
# Why this matters: the third invariant in this file's docstring is "exactly
# one SSE fanout per replayed op." If `with_idempotency`'s cache-replay path
# ever stops short-circuiting the broadcast (e.g. a refactor removes the
# `result.cached` skip), this helper is the only thing that catches it —
# byte-equal HTTP body and durable-row count would both still pass.
function _capture_sse_during(f::Function, kind_filter::String)
    pending = Channel{String}(64)
    sub = (pending = pending,)
    lock(HimalayaUI.SSE_LOCK) do
        push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
    end
    try
        f()
        # Allow the post-commit broadcast queue to flush. `apply_event!` enqueues
        # the broadcast; it fires from the post-commit hook on the same thread.
        sleep(0.3)
    finally
        lock(HimalayaUI.SSE_LOCK) do
            filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
        end
        close(pending)
    end
    frames = String[]
    for frame in pending
        # Skip heartbeats (": heartbeat" lines) and frames of other kinds.
        startswith(frame, ":") && continue
        occursin("\"kind\":\"$kind_filter\"", frame) && push!(frames, frame)
    end
    frames
end

@testset "Idempotency replay invariant — same client_op_id → 1 row, 1 SSE, identical body" begin

    @testset "POST /api/exposures/:id/peaks (peak_added)" begin
        mktempdir() do tmp
            analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
            mkpath(analysis_dir)
            cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
               joinpath(analysis_dir, "example_tot.dat"))
            db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
            exp_id = HimalayaUI.init_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
            s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
            e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
            HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

            with_test_server(db) do port, base
                op_id = "replay-test-peak-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice",
                           "X-Client-Id"  => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:q => 0.99))

                pre_count = _count_actions(db, "peak_added")
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("peak_added") do
                    r1 = HTTP.post("$base/api/exposures/$e_id/peaks";
                        body = body_json, headers = headers)
                    # Replay — same op_id, must produce identical body and zero new rows.
                    r2 = HTTP.post("$base/api/exposures/$e_id/peaks";
                        body = body_json, headers = headers)
                end
                @test r1.status == 201
                @test r2.status == 201
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(db, "peak_added")
                @test post_count - pre_count == 1  # exactly one durable row
                @test length(frames) == 1          # exactly one SSE fanout
            end
        end
    end

    @testset "POST /api/groups/:id/members (index_confirmed)" begin
        mktempdir() do tmp
            analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
            mkpath(analysis_dir)
            cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
               joinpath(analysis_dir, "example_tot.dat"))
            db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
            exp_id = HimalayaUI.init_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
            s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
            e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
            HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

            # Get an index id and the auto group id.
            ix_id = Int(Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM indices WHERE exposure_id = ? LIMIT 1",
                [e_id]))[1].id)
            grp_rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM index_groups WHERE exposure_id = ? AND kind = 'custom'",
                [e_id]))
            isempty(grp_rows) && (HimalayaUI.ensure_custom_group!(db, e_id);
                grp_rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id FROM index_groups WHERE exposure_id = ? AND kind = 'custom'",
                    [e_id])))
            grp_id = Int(grp_rows[1].id)

            with_test_server(db) do port, base
                op_id = "replay-test-grp-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice",
                           "X-Client-Id"  => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:index_id => ix_id))

                pre_count = _count_actions(db, "index_confirmed")
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("index_confirmed") do
                    r1 = HTTP.post("$base/api/groups/$grp_id/members";
                        body = body_json, headers = headers)
                    r2 = HTTP.post("$base/api/groups/$grp_id/members";
                        body = body_json, headers = headers)
                end
                @test r1.status == 200
                @test r2.status == 200
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(db, "index_confirmed")
                @test post_count - pre_count == 1
                @test length(frames) == 1
            end
        end
    end

    @testset "POST /api/samples/:id/messages (post_message)" begin
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
            s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")

            with_test_server(db) do port, base
                op_id = "replay-test-msg-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice",
                           "X-Client-Id"  => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:body => "hello"))

                pre_count = _count_actions(db, "post_message")
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("post_message") do
                    r1 = HTTP.post("$base/api/samples/$s_id/messages";
                        body = body_json, headers = headers)
                    r2 = HTTP.post("$base/api/samples/$s_id/messages";
                        body = body_json, headers = headers)
                end
                @test r1.status == 201
                @test r2.status == 201
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(db, "post_message")
                @test post_count - pre_count == 1
                @test length(frames) == 1
            end
        end
    end

    @testset "PATCH /api/samples/:id (update_sample)" begin
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
            s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")

            with_test_server(db) do port, base
                op_id = "replay-test-upd-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice",
                           "X-Client-Id"  => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:name => "renamed"))

                pre_count = _count_actions(db, "update_sample")
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("update_sample") do
                    r1 = HTTP.patch("$base/api/samples/$s_id";
                        body = body_json, headers = headers)
                    r2 = HTTP.patch("$base/api/samples/$s_id";
                        body = body_json, headers = headers)
                end
                @test r1.status == 200
                @test r2.status == 200
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(db, "update_sample")
                @test post_count - pre_count == 1
                @test length(frames) == 1
            end
        end
    end

    # ── Compare page: comparison_created / comparison_submitted /
    # comparison_deleted SSE-fanout-under-retry rows. Each event kind gets
    # its own row so a regression in `with_idempotency`'s broadcast skip
    # path is caught for each route's enqueue site.

    function _replay_setup_compare(tmp::String)
        analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
        mkpath(analysis_dir)
        cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
           joinpath(analysis_dir, "example_tot.dat"))
        db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
        exp_id = HimalayaUI.init_experiment!(db; path=tmp,
            data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
        s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
        e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
        HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
        (db = db, exposure_id = e_id)
    end

    function _create_body_replay(eid::Int; title::String = "T")
        Dict{Symbol,Any}(:title => title,
            :members => [Dict(:exposure_id => eid, :display_order => 0,
                              :snapshot => Dict(:effective_peaks => [],
                                                :confirmed_index => nothing,
                                                :analysis_inputs_hash => "sha256:zero"))])
    end

    @testset "POST /api/comparisons (comparison_created)" begin
        mktempdir() do tmp
            ctx = _replay_setup_compare(tmp)
            with_test_server(ctx.db) do port, base
                op_id = "replay-cmp-create-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username" => "alice",
                           "X-Client-Id" => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(_create_body_replay(ctx.exposure_id))

                pre_count = _count_actions(ctx.db, "comparison_created")
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("comparison_created") do
                    r1 = HTTP.post("$base/api/comparisons";
                        body = body_json, headers = headers)
                    r2 = HTTP.post("$base/api/comparisons";
                        body = body_json, headers = headers)
                end
                @test r1.status == 201
                @test r2.status == 201
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(ctx.db, "comparison_created")
                @test post_count - pre_count == 1
                @test length(frames) == 1
            end
        end
    end

    @testset "POST /api/comparisons/:id/submit (comparison_submitted)" begin
        mktempdir() do tmp
            ctx = _replay_setup_compare(tmp)
            with_test_server(ctx.db) do port, base
                # Create first.
                r0 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body_replay(ctx.exposure_id; title="orig")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username" => "alice"])
                cmp = JSON3.read(String(r0.body))
                m_id = cmp.members[1].id

                op_id = "replay-cmp-submit-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username" => "alice",
                           "X-Client-Id" => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:title => "renamed",
                    :description => nothing,
                    :expected_content_hash => cmp.content_hash,
                    :members => [Dict(:id => m_id, :exposure_id => ctx.exposure_id,
                                      :display_order => 0,
                                      :snapshot => Dict(:effective_peaks => [],
                                                        :confirmed_index => nothing,
                                                        :analysis_inputs_hash => "sha256:zero"))]))

                pre_count = _count_actions(ctx.db, "comparison_submitted")
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("comparison_submitted") do
                    r1 = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                        body = body_json, headers = headers)
                    r2 = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                        body = body_json, headers = headers)
                end
                @test r1.status == 200
                @test r2.status == 200
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(ctx.db, "comparison_submitted")
                @test post_count - pre_count == 1
                @test length(frames) == 1
            end
        end
    end

    @testset "DELETE /api/comparisons/:id (comparison_deleted)" begin
        mktempdir() do tmp
            ctx = _replay_setup_compare(tmp)
            with_test_server(ctx.db) do port, base
                r0 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body_replay(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username" => "alice"])
                cmp = JSON3.read(String(r0.body))

                op_id = "replay-cmp-del-$(rand(UInt32))"
                headers = ["X-Username" => "alice",
                           "X-Client-Id" => "tab-1",
                           "X-Client-Op-Id" => op_id]

                pre_count = _count_actions(ctx.db, "comparison_deleted")
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("comparison_deleted") do
                    r1 = HTTP.delete("$base/api/comparisons/$(cmp.id)"; headers = headers)
                    r2 = HTTP.delete("$base/api/comparisons/$(cmp.id)"; headers = headers)
                end
                @test r1.status == 200
                @test r2.status == 200
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(ctx.db, "comparison_deleted")
                @test post_count - pre_count == 1
                @test length(frames) == 1
            end
        end
    end

    @testset "POST /api/comparisons/:id/messages (post_message comparison_message)" begin
        mktempdir() do tmp
            ctx = _replay_setup_compare(tmp)
            with_test_server(ctx.db) do port, base
                r0 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body_replay(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username" => "alice"])
                cmp = JSON3.read(String(r0.body))

                op_id = "replay-cmp-msg-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username" => "alice",
                           "X-Client-Id" => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:body => "hello"))

                # Filter post_message but only the comparison_message kind. We
                # rely on per-frame `entity_type` to discriminate from the
                # sample_messages flavor.
                pre_count = first(Tables.rowtable(DBInterface.execute(ctx.db,
                    """SELECT COUNT(*) AS c FROM user_actions
                       WHERE action = 'post_message' AND entity_type = 'comparison_message'"""))).c
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("post_message") do
                    r1 = HTTP.post("$base/api/comparisons/$(cmp.id)/messages";
                        body = body_json, headers = headers)
                    r2 = HTTP.post("$base/api/comparisons/$(cmp.id)/messages";
                        body = body_json, headers = headers)
                end
                @test r1.status == 201
                @test r2.status == 201
                @test String(r2.body) == String(r1.body)
                post_count = first(Tables.rowtable(DBInterface.execute(ctx.db,
                    """SELECT COUNT(*) AS c FROM user_actions
                       WHERE action = 'post_message' AND entity_type = 'comparison_message'"""))).c
                @test post_count - pre_count == 1
                # Filter to only comparison_message frames (same kind name as
                # sample_messages but discriminated by entity_type in the JSON).
                comp_frames = filter(f -> occursin("\"entity_type\":\"comparison_message\"", f),
                                     frames)
                @test length(comp_frames) == 1
            end
        end
    end

    @testset "POST /api/comparisons/:id/pin (comparison_pinned)" begin
        mktempdir() do tmp
            ctx = _replay_setup_compare(tmp)
            with_test_server(ctx.db) do port, base
                # Pre-mint a comparison via the route so the pin target exists.
                r0 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body_replay(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username" => "alice"])
                cmp = JSON3.read(String(r0.body))

                op_id = "replay-cmp-pin-$(rand(UInt32))"
                headers = ["X-Username" => "alice",
                           "X-Client-Id" => "tab-1",
                           "X-Client-Op-Id" => op_id]

                pre_count = _count_actions(ctx.db, "comparison_pinned")
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("comparison_pinned") do
                    r1 = HTTP.post("$base/api/comparisons/$(cmp.id)/pin"; headers = headers)
                    r2 = HTTP.post("$base/api/comparisons/$(cmp.id)/pin"; headers = headers)
                end
                @test r1.status == 200
                @test r2.status == 200
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(ctx.db, "comparison_pinned")
                @test post_count - pre_count == 1
                @test length(frames) == 1
            end
        end
    end

    @testset "DELETE /api/comparisons/:id/pin (comparison_unpinned)" begin
        mktempdir() do tmp
            ctx = _replay_setup_compare(tmp)
            with_test_server(ctx.db) do port, base
                r0 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body_replay(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username" => "alice"])
                cmp = JSON3.read(String(r0.body))
                # Establish a pin to unpin.
                HTTP.post("$base/api/comparisons/$(cmp.id)/pin",
                    ["X-Username" => "alice"])

                op_id = "replay-cmp-unpin-$(rand(UInt32))"
                headers = ["X-Username" => "alice",
                           "X-Client-Id" => "tab-1",
                           "X-Client-Op-Id" => op_id]

                pre_count = _count_actions(ctx.db, "comparison_unpinned")
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("comparison_unpinned") do
                    r1 = HTTP.delete("$base/api/comparisons/$(cmp.id)/pin"; headers = headers)
                    r2 = HTTP.delete("$base/api/comparisons/$(cmp.id)/pin"; headers = headers)
                end
                @test r1.status == 200
                @test r2.status == 200
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(ctx.db, "comparison_unpinned")
                @test post_count - pre_count == 1
                @test length(frames) == 1
            end
        end
    end

    @testset "POST /submit 409 retry: status >= 400 NOT cached → conflict re-evaluates" begin
        mktempdir() do tmp
            ctx = _replay_setup_compare(tmp)
            with_test_server(ctx.db) do port, base
                r0 = HTTP.post("$base/api/comparisons";
                    body = JSON3.write(_create_body_replay(ctx.exposure_id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username" => "alice"])
                cmp = JSON3.read(String(r0.body))
                m_id = cmp.members[1].id

                op_id = "conflict-replay-$(rand(UInt32))"
                stale_body = JSON3.write(Dict(:title => "x",
                    :description => nothing,
                    :expected_content_hash => "sha256:stale",
                    :members => [Dict(:id => m_id, :exposure_id => ctx.exposure_id,
                                      :display_order => 0,
                                      :snapshot => Dict(:effective_peaks => [],
                                                        :confirmed_index => nothing,
                                                        :analysis_inputs_hash => "sha256:zero"))]))
                hdrs = ["Content-Type" => "application/json",
                        "X-Username" => "alice",
                        "X-Client-Op-Id" => op_id]

                # Both calls return 409 (re-evaluated). 4xx not cached →
                # the second call walks the same body, which still finds a
                # stale expected_content_hash and 409s again. NO durable
                # `comparison_submitted` row is written for either.
                pre_submitted = _count_actions(ctx.db, "comparison_submitted")
                r1 = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                    body = stale_body, headers = hdrs, status_exception = false)
                r2 = HTTP.post("$base/api/comparisons/$(cmp.id)/submit";
                    body = stale_body, headers = hdrs, status_exception = false)
                @test r1.status == 409
                @test r2.status == 409
                # Both bodies present a fresh current_state — they are NOT
                # byte-identical cached replays (the conflict was re-evaluated).
                # current_hash is the same because the underlying state didn't
                # change between calls; but the retry actually ran the route.
                post_submitted = _count_actions(ctx.db, "comparison_submitted")
                @test post_submitted == pre_submitted
            end
        end
    end
end
