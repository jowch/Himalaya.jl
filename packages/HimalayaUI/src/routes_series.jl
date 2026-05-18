using HTTP, JSON3, DBInterface, Tables, Oxygen, SQLite

# ─────────────────────────────────────────────────────────────────────────────
# I2.2 (#165) — REST routes for the series model. Mirrors
# `routes_comparisons.jl` (corpus-wide — there is no `/api/experiments/{eid}/series`).
#
# Every mutating route wraps in `with_idempotency(db, req)` and uses
# `apply_event!(InTransaction(), …)`. Body-shape validation runs BEFORE
# `with_idempotency` so a malformed payload returns an uncached 400.
#
# No route carries an `is_author` / 403 gate (architecture decision 3).
# Existence (404) and optimistic-concurrency (409) checks remain.
#
# `_json_error` and `_view_fields_error` are reused from `routes_comparisons.jl`
# (same module). Phase 3 (#175 / I3.6) relocates them when it deletes that file.
# ─────────────────────────────────────────────────────────────────────────────

function register_series_routes!()
    # ── Listing ─────────────────────────────────────────────────────────────

    @get "/api/series" function(req::HTTP.Request)
        db = current_db()
        rows = series_listing(db)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(rows))
    end

    # ── Detail ──────────────────────────────────────────────────────────────

    @get "/api/series/{id}" function(req::HTTP.Request, id::Int)
        db = current_db()
        out = fetch_series_with_plate(db, id)
        out === nothing && return _json_error(404, "series not found")
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
    end

    @get "/api/series/{id}/forks" function(req::HTTP.Request, id::Int)
        db = current_db()
        rows = forks_of_series(db, id)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(rows))
    end

    # ── Create (draft) ──────────────────────────────────────────────────────

    @post "/api/series" function(req::HTTP.Request)
        db   = current_db()
        body = json(req)

        # Validate the body shape BEFORE with_idempotency so a malformed
        # payload returns an uncached 400 (validation toast), not a cached 500.
        if !haskey(body, :title)
            return _json_error(400, "missing required field: title")
        end
        if haskey(body, :samples) && body.samples !== nothing
            if !(body.samples isa AbstractVector)
                return _json_error(400, "samples must be an array")
            end
            # `sample_id` is required per entry — `_series_sample_payload`
            # indexes it unconditionally, and a recipe row with no target is
            # unrenderable. Reject here so a bad entry is an uncached 400, not
            # a 500 cached inside with_idempotency.
            for m in body.samples
                if !(m isa AbstractDict || m isa JSON3.Object) ||
                        !haskey(m, :sample_id) || m.sample_id === nothing
                    return _json_error(400, "each samples entry requires sample_id")
                end
            end
        end
        verr = _view_fields_error(body)
        verr === nothing || return verr

        title = String(body.title)
        description = haskey(body, :description) && body.description !== nothing ?
                      String(body.description) : nothing
        forked_from_id = haskey(body, :forked_from_id) && body.forked_from_id !== nothing ?
                         Int(body.forked_from_id) : nothing
        forked_at_hash = haskey(body, :forked_at_hash) && body.forked_at_hash !== nothing ?
                         String(body.forked_at_hash) : nothing
        ordering_variable = haskey(body, :ordering_variable) && body.ordering_variable !== nothing ?
                            String(body.ordering_variable) : nothing
        order_rule = haskey(body, :order_rule) && body.order_rule !== nothing ?
                     String(body.order_rule) : nothing
        view_grouping_mode = haskey(body, :view_grouping_mode) && body.view_grouping_mode !== nothing ?
                             String(body.view_grouping_mode) : nothing
        # Boolean view fields: `haskey` only (no `!== nothing` guard) — `false`
        # is a valid value. An explicit JSON `null` passes through as `nothing`,
        # which `_view_fields_error` permits as "reset to the per-tab default".
        view_show_peak_ticks  = haskey(body, :view_show_peak_ticks)  ? body.view_show_peak_ticks  : nothing
        view_show_peak_labels = haskey(body, :view_show_peak_labels) ? body.view_show_peak_labels : nothing
        samples_in = haskey(body, :samples) && body.samples !== nothing ? body.samples : ()

        return with_idempotency(db, req) do
            # Mint the AUTOINCREMENT id with a NULL-only placeholder row. The
            # `series_created` dispatcher (#166) upserts at this id and sets
            # `state='draft'`; until it lands the placeholder keeps the schema
            # default `state='committed'` — accepted (the row is degenerate).
            res = DBInterface.execute(db, "INSERT INTO series DEFAULT VALUES")
            new_id = Int(DBInterface.lastrowid(res))

            samples_payload = [_series_sample_payload(m) for m in samples_in]
            payload = Dict{Symbol, Any}(
                :title                 => title,
                :description           => description,
                :forked_from_id        => forked_from_id,
                :forked_at_hash        => forked_at_hash,
                :ordering_variable     => ordering_variable,
                :order_rule            => order_rule,
                :view_grouping_mode    => view_grouping_mode,
                :view_show_peak_ticks  => view_show_peak_ticks,
                :view_show_peak_labels => view_show_peak_labels,
                :samples               => samples_payload,
            )
            result = apply_event!(InTransaction(), db, req;
                kind        = "series_created",
                entity_type = "series",
                entity_id   = new_id,
                payload     = payload)

            out = fetch_series_with_plate(db, new_id)
            out === nothing && error(
                "post-write fetch_series_with_plate returned nothing for id=$(new_id)")
            # series_created carries NO post_state envelope (master-plan §5.2);
            # foreign tabs reconcile via replay-as-rerun (#166).
            _enqueue_broadcast_from_result!(result, "series_created",
                                            "series", new_id)
            HTTP.Response(201, ["Content-Type" => "application/json"], JSON3.write(out))
        end
    end
end
