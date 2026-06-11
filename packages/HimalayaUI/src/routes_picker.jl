using HTTP, JSON3, DBInterface, Tables, Oxygen

# ─────────────────────────────────────────────────────────────────────────────
# Phase 5 Task 5.2 — picker support routes.
#
# Two read-only GETs feed the frontend ComparisonPicker modal:
#
#   GET /api/users/:id/recently-picked-exposures?limit=20
#       The user's most-recently-picked exposures, derived from
#       `comparison_members.created_by/created_at` via
#       `recently_used_exposures`. Surfaces the per-user history that the
#       picker uses to populate its "Recently used" section.
#
#   GET /api/experiments/:eid/sample-tags
#       Distinct (key, value) pairs across every sample in the experiment,
#       used to populate the picker's tag-filter dropdown. CRITICAL: distinct
#       collapses on the (key, value) pair, NOT on value alone — two tags
#       with the same value but different keys must surface as separate
#       options. The picker treats each (key, value) as an independent
#       filter chip per the spec.
#
#   GET /api/sample-tags
#       Corpus-wide sibling of the route above: distinct (key, value) pairs
#       across every sample in the whole database, with no experiment
#       scoping. Feeds the redesign's corpus-level series scoping step
#       (issue map I0.2). Same (key, value) distinct contract.
#
# All three routes are read-only — no `with_idempotency`, no SSE, no event-log
# row. Callers are not expected to mutate state through these endpoints.
# ─────────────────────────────────────────────────────────────────────────────

function register_picker_routes!()
    @get "/api/users/{id}/recently-picked-exposures" function(req::HTTP.Request, id::Int)
        db = current_db()
        # Validate the user exists; return 404 (not an empty list) for an
        # unknown id so the frontend can distinguish "nothing yet" from
        # "wrong url". Mirrors the `/api/users/:username/actions` 404 path.
        urows = Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM users WHERE id = ?", [id]))
        isempty(urows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "user not found")))

        # Limit query param: positive int, default 20, clamped at 100 so a
        # malicious caller can't ask for an unbounded list.
        limit = 20
        params = HTTP.queryparams(HTTP.URI(req.target))
        if haskey(params, "limit")
            try
                lim = parse(Int, params["limit"])
                limit = clamp(lim, 1, 100)
            catch
                # Bad limit → fall back to default; invalid input shouldn't 500.
            end
        end

        ids = recently_used_exposures(db, id; limit = limit)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(ids))
    end

    @get "/api/experiments/{eid}/sample-tags" function(req::HTTP.Request, eid::Int)
        db = current_db()
        # Per-(key, value) sample count joined through samples → experiment.
        # No need to validate the experiment id — an unknown id naturally
        # produces an empty list, which is the right answer for the picker.
        rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT t.key, t.value, COUNT(DISTINCT t.sample_id) AS count
               FROM sample_tags t
               JOIN samples s ON s.id = t.sample_id
               WHERE s.experiment_id = ?
               GROUP BY t.key, t.value
               ORDER BY t.key, t.value""", [eid]))

        out = [Dict(:key => String(r.key), :value => String(r.value), :count => Int(r.count)) for r in rows]
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
    end

    @get "/api/sample-tags" function(req::HTTP.Request)
        db = current_db()
        # Per-(key, value) sample count across the whole corpus, with no
        # experiment filter. COUNT(DISTINCT sample_id) gives the number of
        # samples carrying each (key, value) pair — used by the Manage-tags
        # modal to rank suggestions by frequency. proposeOrdering (scoping)
        # ranks by distinct-value count instead and ignores this field.
        rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT key, value, COUNT(DISTINCT sample_id) AS count
               FROM sample_tags
               GROUP BY key, value
               ORDER BY key, value"""))

        out = [Dict(:key => String(r.key), :value => String(r.value), :count => Int(r.count)) for r in rows]
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
    end

    @get "/api/experiments/{eid}/picker-samples" function(req::HTTP.Request, eid::Int)
        db = current_db()
        rows = picker_samples(db, eid)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(rows))
    end

    @get "/api/picker-samples" function(req::HTTP.Request)
        # Corpus-wide sibling of the experiment-scoped picker-samples route.
        # Read-only — no with_idempotency, no SSE, no event-log row, matching
        # the other picker routes.
        db = current_db()
        rows = picker_samples(db)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(rows))
    end
end
