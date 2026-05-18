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
end
