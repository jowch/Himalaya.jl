using HTTP, JSON3, DBInterface, Tables, Oxygen, SQLite

# ─────────────────────────────────────────────────────────────────────────────
# Phase 2 Task 2.2 — REST routes for the Compare page.
#
# Every mutating route wraps in `with_idempotency(db, req)` and uses
# `apply_event!(InTransaction(), …)` so the durable event row, the view-table
# update, and the idempotency cache row commit (or roll back) atomically.
# SSE broadcasts fire post-commit via the queued-broadcast helper.
#
# Author-only gates use `is_author(db, comparison_id, user_id)` which returns
# false for missing comparisons, NULL `created_by` (orphan author), and any
# user mismatch — see comparisons.jl docstring.
# ─────────────────────────────────────────────────────────────────────────────

"""
    _comparison_member_payload(db, m_in)

Normalize a member entry from the request body into the payload shape the
dispatcher expects. Fills `snapshot` from `compute_member_snapshot(db, …)`
when the client omitted it, AND when an `exposure_id` is present (orphan
members have no live snapshot to compute). The dispatcher errors on a
member with a missing snapshot, so this normalization is load-bearing for
the create-without-client-snapshot path.
"""
function _comparison_member_payload(db::SQLite.DB, m_in)
    eid = haskey(m_in, :exposure_id) ? m_in[:exposure_id] : nothing
    snap_in = haskey(m_in, :snapshot) ? m_in[:snapshot] : nothing
    snap = if snap_in !== nothing
        snap_in
    elseif eid !== nothing
        compute_member_snapshot(db, Int(eid))
    else
        # Orphan member with no snapshot — fall back to a minimal one so the
        # NOT NULL CHECK in comparison_members.snapshot is satisfied.
        Dict{Symbol, Any}(
            :effective_peaks      => Any[],
            :confirmed_index      => nothing,
            :analysis_inputs_hash => nothing,
        )
    end

    out = Dict{Symbol, Any}(
        :id             => haskey(m_in, :id) ? m_in[:id] : nothing,
        :exposure_id    => eid,
        :display_order  => haskey(m_in, :display_order) ? Int(m_in[:display_order]) : 0,
        :band_height    => haskey(m_in, :band_height)   ? Float64(m_in[:band_height]) : 1.0,
        :y_offset       => haskey(m_in, :y_offset)      ? Float64(m_in[:y_offset])    : 0.0,
        :normalization  => haskey(m_in, :normalization) ? String(m_in[:normalization]) : "none",
        :color_override => haskey(m_in, :color_override) ? m_in[:color_override] : nothing,
        :label_override => haskey(m_in, :label_override) ? m_in[:label_override] : nothing,
        :q_window_min   => haskey(m_in, :q_window_min)   ? m_in[:q_window_min]   : nothing,
        :q_window_max   => haskey(m_in, :q_window_max)   ? m_in[:q_window_max]   : nothing,
        :peak_display   => haskey(m_in, :peak_display)   ? m_in[:peak_display]   : nothing,
        :snapshot       => snap,
    )
    out
end

"""
    _json_error(status, msg; extra...)

Single-shot JSON error response with `Content-Type: application/json`.
"""
function _json_error(status::Int, msg::AbstractString; extra...)
    body = Dict{Symbol, Any}(:error => msg)
    for (k, v) in extra
        body[k] = v
    end
    HTTP.Response(status, ["Content-Type" => "application/json"],
                  JSON3.write(body))
end

function register_comparisons_routes!()
    # ── Listing routes ──────────────────────────────────────────────────────

    @get "/api/experiments/{eid}/comparisons" function(req::HTTP.Request, eid::Int)
        db = current_db()
        rows = comparisons_for_experiment(db, eid)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(rows))
    end

    @get "/api/comparisons" function(req::HTTP.Request)
        db = current_db()
        rows = comparisons_listing(db)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(rows))
    end

    # ── Create + initial submit ─────────────────────────────────────────────

    @post "/api/comparisons" function(req::HTTP.Request)
        db   = current_db()
        body = json(req)

        # Validate body shape BEFORE entering with_idempotency so a malformed
        # payload returns 400 (validation toast on the frontend) rather than
        # 500 (infrastructure banner). 400s are not cached, so retries with
        # the same client_op_id re-evaluate.
        if !haskey(body, :title)
            return _json_error(400, "missing required field: title")
        end
        if !haskey(body, :members) || !(body.members isa AbstractVector) ||
                isempty(body.members)
            return _json_error(400, "members must be a non-empty array")
        end

        title = String(body.title)
        description = haskey(body, :description) && body.description !== nothing ?
                      String(body.description) : nothing
        forked_from_id = haskey(body, :forked_from_id) && body.forked_from_id !== nothing ?
                         Int(body.forked_from_id) : nothing
        forked_at_hash = haskey(body, :forked_at_hash) && body.forked_at_hash !== nothing ?
                         String(body.forked_at_hash) : nothing

        return with_idempotency(db, req) do
            # Mint the AUTOINCREMENT id with a NULL-only placeholder row.
            # The dispatcher's `comparison_created` branch upserts at this id
            # via plain `COALESCE(col, ?)` — see
            # _update_view_for_comparison_created!. Pre-#67 this INSERT seeded
            # `''` for four NOT NULL columns; #67 relaxed those to nullable
            # so `DEFAULT VALUES` Just Works and no NULLIF wrapper is needed
            # at the dispatcher.
            res = DBInterface.execute(db,
                "INSERT INTO comparisons DEFAULT VALUES")
            new_id = Int(DBInterface.lastrowid(res))

            members_payload = [_comparison_member_payload(db, m) for m in body.members]
            # Strip ids from create-mode members — IDs are minted by the
            # dispatcher on INSERT (the `comparison_created` branch ignores
            # any payload ids per its docstring).
            for m in members_payload
                m[:id] = nothing
            end

            payload = Dict{Symbol, Any}(
                :title          => title,
                :description    => description,
                :forked_from_id => forked_from_id,
                :forked_at_hash => forked_at_hash,
                :members        => members_payload,
            )
            result = apply_event!(InTransaction(), db, req;
                kind        = "comparison_created",
                entity_type = "comparison",
                entity_id   = new_id,
                payload     = payload)
            _enqueue_broadcast_from_result!(result, "comparison_created",
                                            "comparison", new_id)

            out = fetch_comparison_with_members(db, new_id)
            HTTP.Response(201, ["Content-Type" => "application/json"], JSON3.write(out))
        end
    end

    # ── Detail + lineage ────────────────────────────────────────────────────

    @get "/api/comparisons/{id}" function(req::HTTP.Request, id::Int)
        db = current_db()
        out = fetch_comparison_with_members(db, id)
        out === nothing && return _json_error(404, "comparison not found")
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
    end

    @get "/api/comparisons/{id}/forks" function(req::HTTP.Request, id::Int)
        db = current_db()
        rows = forks_of_comparison(db, id)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(rows))
    end

    # ── Submit (author-only, conflict-aware) ────────────────────────────────

    @post "/api/comparisons/{id}/submit" function(req::HTTP.Request, id::Int)
        db = current_db()
        body = json(req)

        if !haskey(body, :title)
            return _json_error(400, "missing required field: title")
        end
        if !haskey(body, :members) || !(body.members isa AbstractVector)
            return _json_error(400, "members must be an array")
        end

        title = String(body.title)
        description = haskey(body, :description) && body.description !== nothing ?
                      String(body.description) : nothing
        expected_hash = haskey(body, :expected_content_hash) &&
                        body.expected_content_hash !== nothing ?
                        String(body.expected_content_hash) : nothing

        return with_idempotency(db, req) do
            # Existence check before the author gate: HTTP semantics require
            # 404 to fire BEFORE 403 — a request for a resource that doesn't
            # exist can't be "forbidden" because there's nothing to forbid.
            # `current_content_hash === nothing` is the existence probe (it
            # returns nothing for missing comparisons).
            if current_content_hash(db, id) === nothing
                return _json_error(404, "comparison not found")
            end
            # Author gate. 403s are not cached → a retry under the same
            # client_op_id re-evaluates (allows the conflict to resolve if
            # ownership changes via fork or admin action).
            user_id = get_user_id_for_request(db, req)
            if !is_author(db, id, user_id)
                # The orphan-author case (`created_by IS NULL`) is covered
                # by `is_author` returning false; the comparison exists but
                # has no author.
                return _json_error(403, "only the author can submit")
            end

            # Conflict check: stored hash MUST match the client's
            # expected_content_hash, otherwise return 409 with current_state.
            current_hash = current_content_hash(db, id)
            if expected_hash !== nothing && current_hash !== expected_hash
                current_state = fetch_comparison_with_members(db, id)
                return HTTP.Response(409, ["Content-Type" => "application/json"],
                    JSON3.write(Dict(
                        :error         => "conflict",
                        :current_hash  => current_hash,
                        :current_state => current_state,
                    )))
            end

            members_payload = [_comparison_member_payload(db, m) for m in body.members]
            payload = Dict{Symbol, Any}(
                :title       => title,
                :description => description,
                :members     => members_payload,
            )
            result = apply_event!(InTransaction(), db, req;
                kind        = "comparison_submitted",
                entity_type = "comparison",
                entity_id   = id,
                payload     = payload)
            _enqueue_broadcast_from_result!(result, "comparison_submitted",
                                            "comparison", id)

            out = fetch_comparison_with_members(db, id)
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
        end
    end

    # ── Delete (author-only) ────────────────────────────────────────────────

    @delete "/api/comparisons/{id}" function(req::HTTP.Request, id::Int)
        db = current_db()
        return with_idempotency(db, req) do
            # Existence check before the author gate (see POST /submit for
            # rationale). `current_content_hash` is the existence probe.
            if current_content_hash(db, id) === nothing
                return _json_error(404, "comparison not found")
            end
            user_id = get_user_id_for_request(db, req)
            if !is_author(db, id, user_id)
                return _json_error(403, "only the author can delete")
            end

            result = apply_event!(InTransaction(), db, req;
                kind        = "comparison_deleted",
                entity_type = "comparison",
                entity_id   = id,
                payload     = Dict{Symbol, Any}(:id => id))
            _enqueue_broadcast_from_result!(result, "comparison_deleted",
                                            "comparison", id)

            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:id => id, :deleted => true,
                                 :event_id => result.event_id)))
        end
    end

    # ── Chat thread ─────────────────────────────────────────────────────────

    @get "/api/comparisons/{id}/messages" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db, """
            SELECT m.id, m.comparison_id, m.author_id,
                   u.username AS author,
                   m.body, m.created_at
            FROM comparison_messages m
            LEFT JOIN users u ON u.id = m.author_id
            WHERE m.comparison_id = ?
            ORDER BY m.created_at ASC, m.id ASC
        """, [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(rows_to_json(rows)))
    end

    @post "/api/comparisons/{id}/messages" function(req::HTTP.Request, id::Int)
        db       = current_db()
        username = get_username(req)
        if username === nothing
            return _json_error(401, "X-Username header required")
        end

        body = json(req)
        text = haskey(body, :body) ? strip(String(body.body)) : ""
        if isempty(text)
            return _json_error(400, "message body required")
        end

        return with_idempotency(db, req) do
            # Mirror routes_messages.jl: the route writes the message row
            # directly (the dispatcher is a no-op for `post_message`).
            author_id = get_or_create_user!(db, username)
            res = DBInterface.execute(db,
                "INSERT INTO comparison_messages (comparison_id, author_id, body) VALUES (?, ?, ?)",
                [id, author_id, text])
            msg_id = Int(DBInterface.lastrowid(res))

            row = Tables.rowtable(DBInterface.execute(db, """
                SELECT m.id, m.comparison_id, m.author_id,
                       u.username AS author,
                       m.body, m.created_at
                FROM comparison_messages m
                LEFT JOIN users u ON u.id = m.author_id
                WHERE m.id = ?
            """, [msg_id]))[1]
            msg_json = row_to_json(row)

            # `entity_type='comparison_message'` differentiates from the
            # `sample_message` path in `routes_messages.jl`. The dispatcher
            # is a no-op for `post_message` (kind sees both entity types as
            # log-only).
            result = apply_event!(InTransaction(), db, req;
                kind        = "post_message",
                entity_type = "comparison_message",
                entity_id   = msg_id,
                payload     = msg_json)
            _enqueue_broadcast_from_result!(result, "post_message",
                                            "comparison_message", msg_id)

            HTTP.Response(201, ["Content-Type" => "application/json"],
                JSON3.write(msg_json))
        end
    end

    # ── Pins (Phase 13 + Phase 4 follow-up) ────────────────────────────────
    #
    # Per-user pinned comparisons. Pin/unpin route through the event log so:
    #   1. Cross-tab fanout works via SSE (a pin in tab A propagates to tab
    #      B before cache stales; previously tabs only saw their own pins).
    #   2. Durable history is recorded for audit / per-user activity views.
    #   3. The route surface looks like every other mutating route — wrapped
    #      in `with_idempotency`, calls `apply_event!(InTransaction(), …)`,
    #      enqueues a post-commit broadcast.
    #
    # Event entity model: `entity_type='user'`, `entity_id=user_id`. The
    # event "happens to" the user (their pin set changes); the affected
    # comparison rides in the payload. This matches the typical durable-log
    # convention where the entity is the noun whose state mutates.
    #
    # The view-table write (INSERT/DELETE on `comparison_pins`) moves into
    # the dispatcher's `comparison_pinned`/`comparison_unpinned` branches —
    # the dispatcher derives user_id from `user_actions WHERE id = event_id`
    # so the route doesn't need to forward it through the payload.
    #
    # GET /api/users/me/comparison-pins stays a plain read; no event involved.

    @get "/api/users/me/comparison-pins" function(req::HTTP.Request)
        db = current_db()
        user_id = get_user_id_for_request(db, req)
        if user_id === nothing
            return _json_error(401, "X-Username header required")
        end
        # Most-recently-pinned first.
        rows = Tables.rowtable(DBInterface.execute(db, """
            SELECT comparison_id FROM comparison_pins
            WHERE user_id = ?
            ORDER BY pinned_at DESC, comparison_id DESC
        """, [user_id]))
        ids = Int[Int(r.comparison_id) for r in rows]
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(ids))
    end

    @post "/api/comparisons/{id}/pin" function(req::HTTP.Request, id::Int)
        db = current_db()
        user_id = get_user_id_for_request(db, req)
        if user_id === nothing
            return _json_error(401, "X-Username header required")
        end
        return with_idempotency(db, req) do
            # Existence check inside with_idempotency for convention parity
            # with submit/delete. 4xx responses aren't cached by the wrapper,
            # so a future pin attempt after re-creation at a new id still
            # gets a fresh existence probe.
            if current_content_hash(db, id) === nothing
                return _json_error(404, "comparison not found")
            end
            payload = Dict{Symbol, Any}(:comparison_id => id)
            result = apply_event!(InTransaction(), db, req;
                kind        = "comparison_pinned",
                entity_type = "user",
                entity_id   = user_id,
                payload     = payload)
            _enqueue_broadcast_from_result!(result, "comparison_pinned",
                                            "user", user_id)

            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:comparison_id => id, :pinned => true)))
        end
    end

    @delete "/api/comparisons/{id}/pin" function(req::HTTP.Request, id::Int)
        db = current_db()
        user_id = get_user_id_for_request(db, req)
        if user_id === nothing
            return _json_error(401, "X-Username header required")
        end
        return with_idempotency(db, req) do
            # Idempotent at the SQL layer — unpinning a never-pinned
            # comparison is a no-op DELETE, so the frontend doesn't need to
            # distinguish "was pinned" from "wasn't pinned" in a single
            # round trip. The event itself is still recorded so the
            # cross-tab broadcast fires either way.
            payload = Dict{Symbol, Any}(:comparison_id => id)
            result = apply_event!(InTransaction(), db, req;
                kind        = "comparison_unpinned",
                entity_type = "user",
                entity_id   = user_id,
                payload     = payload)
            _enqueue_broadcast_from_result!(result, "comparison_unpinned",
                                            "user", user_id)

            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:comparison_id => id, :pinned => false)))
        end
    end
end
