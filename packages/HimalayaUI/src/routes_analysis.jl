using HTTP, JSON3, DBInterface, Tables, Oxygen, SQLite
using Himalaya

"""
    ensure_custom_group!(db, exposure_id) -> (group_id, created)

Returns the id of the custom group for this exposure and whether it was just
created. If a custom group already exists, returns (existing_id, false).
Otherwise clones the auto group's members into a new custom group (active),
demotes the auto group to inactive, returns (new_id, true).

Errors if no auto group exists for the exposure.
"""
function ensure_custom_group!(db::SQLite.DB, exposure_id::Int)
    existing = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM index_groups
         WHERE exposure_id = ? AND kind = 'custom'", [exposure_id]))
    isempty(existing) || return (Int(existing[1].id), false)

    auto_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM index_groups
         WHERE exposure_id = ? AND kind = 'auto'", [exposure_id]))
    isempty(auto_rows) && error("no auto group for exposure $exposure_id")
    auto_id = Int(auto_rows[1].id)

    res = DBInterface.execute(db,
        "INSERT INTO index_groups (exposure_id, kind, active)
         VALUES (?, 'custom', 1)", [exposure_id])
    custom_id = Int(DBInterface.lastrowid(res))

    DBInterface.execute(db,
        "INSERT INTO index_group_members (group_id, index_id)
         SELECT ?, index_id FROM index_group_members WHERE group_id = ?",
        [custom_id, auto_id])

    DBInterface.execute(db,
        "UPDATE index_groups SET active = 0 WHERE id = ?", [auto_id])

    (custom_id, true)
end

function _group_with_members(db::SQLite.DB, group_id::Int)
    g = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM index_groups WHERE id = ?", [group_id]))[1]
    members = Tables.rowtable(DBInterface.execute(db,
        "SELECT index_id FROM index_group_members
         WHERE group_id = ? ORDER BY index_id", [group_id]))
    d = row_to_json(g; bool_keys = (:active,))
    d[:members] = [Int(m.index_id) for m in members]
    d
end

"""
    predicted_q_for_phase(phase_name, basis) -> Vector{Float64}

Return predicted q positions (basis × normalized phase ratios) for a phase
named by its string form (e.g., "Pn3m"). Returns empty vector if the phase
is unknown to Himalaya.
"""
function predicted_q_for_phase(phase_name::AbstractString, basis::Float64)::Vector{Float64}
    # Strip module prefix if present (e.g. "Himalaya.Pn3m" -> "Pn3m")
    bare = last(split(phase_name, '.'))
    P = try
        getfield(Himalaya, Symbol(bare))
    catch
        return Float64[]
    end
    (P isa Type && P <: Himalaya.Phase) || return Float64[]
    ratios = Himalaya.phaseratios(P; normalize=true)
    [basis * r for r in ratios]
end

"""
    _ngc_for_phase(phase_name, lattice_d) -> Union{Float64, Nothing}

Compute phase-averaged Gaussian curvature κ for the three bicontinuous
cubic phases (Ia3d, Pn3m, Im3m). Result has units `1/length²` matching
the unit of `lattice_d` — i.e. if the experiment stores q in `A-1`,
lattice_d is in Å and κ comes back in `Å⁻²`; if q is in `nm-1`, κ is
in `nm⁻²`. Returns `nothing` for other phases or when `lattice_d` is
missing/non-positive.

Formula: `κ = -2π·(χ/A₀)/a²`. The legacy `Himalaya.ngc` kernel hard-codes
a `(10/a)²` factor that assumed `a` in Å and produced `nm⁻²` output —
fine for the matlab UI it was ported from, but it breaks dimensional
honesty when the experiment isn't in Å. The unit-naked formula here lets
the frontend render κ alongside its experiment-native lattice unit
(`Å⁻²` or `nm⁻²`) without backend-side conversion.
"""
function _ngc_for_phase(phase_name::AbstractString, lattice_d)::Union{Float64, Nothing}
    (lattice_d === nothing || ismissing(lattice_d)) && return nothing
    a = Float64(lattice_d)
    a > 0 || return nothing
    bare = last(split(phase_name, '.'))
    χ, A₀ = if bare == "Ia3d"
        (-8, 3.091)
    elseif bare == "Pn3m"
        (-2, 1.919)
    elseif bare == "Im3m"
        (-4, 2.345)
    else
        return nothing
    end
    -2π * (χ / A₀) / a^2
end

function register_analysis_routes!()
    @get "/api/exposures/{id}/indices" function(req::HTTP.Request, id::Int)
        db = current_db()
        indices = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM indices WHERE exposure_id = ? ORDER BY score DESC", [id]))
        out = map(indices) do ix
            peak_rows = Tables.rowtable(DBInterface.execute(db,
                """SELECT ip.peak_id, ip.ratio_position, ip.residual,
                          COALESCE(ap.q, pc.q) AS q_observed
                   FROM index_peaks ip
                   LEFT JOIN auto_peaks ap     ON ap.id = ip.peak_id AND ip.peak_kind = 'auto'
                   LEFT JOIN peak_curations pc ON pc.id = ip.peak_id AND ip.peak_kind = 'curation'
                   WHERE ip.index_id = ? ORDER BY ip.ratio_position""",
                [Int(ix.id)]))
            predicted = predicted_q_for_phase(String(ix.phase), Float64(ix.basis))
            d = row_to_json(ix)
            d[:peaks]        = rows_to_json(peak_rows)
            d[:predicted_q]  = predicted
            d[:ngc]          = _ngc_for_phase(String(ix.phase), ix.lattice_d)
            d
        end
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(out))
    end

    @get "/api/exposures/{id}/groups" function(req::HTTP.Request, id::Int)
        db     = current_db()
        groups = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM index_groups WHERE exposure_id = ? ORDER BY id", [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write([_group_with_members(db, Int(g.id)) for g in groups]))
    end

    @post "/api/groups/{id}/members" function(req::HTTP.Request, id::Int)
        db = current_db()
        return with_idempotency(db, req) do
            body = json(req)
            index_id = Int(body.index_id)

            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT exposure_id, kind FROM index_groups WHERE id = ?", [id]))
            isempty(rows) && return HTTP.Response(404,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "group not found")))
            exposure_id = Int(rows[1].exposure_id)

            custom_id, _ = ensure_custom_group!(db, exposure_id)

            result = apply_event!(InTransaction(), db, req;
                kind        = "index_confirmed",
                entity_type = "exposure",
                entity_id   = exposure_id,
                payload     = Dict(:group_id => custom_id, :index_id => index_id))
            _enqueue_broadcast_from_result!(result, "index_confirmed", "exposure", exposure_id)

            # Issue #13: include event_id/view_row_id alongside the group
            # body so the response shape matches the spec's queue-migrated
            # contract (event_id, view_row_id, ...).
            body = _group_with_members(db, custom_id)
            body[:event_id]    = result.event_id
            body[:view_row_id] = result.view_row_id
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(body))
        end
    end

    @delete "/api/groups/{id}/members/{index_id}" function(req::HTTP.Request, id::Int, index_id::Int)
        db = current_db()
        return with_idempotency(db, req) do
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT exposure_id FROM index_groups WHERE id = ?", [id]))
            isempty(rows) && return HTTP.Response(404,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "group not found")))
            exposure_id = Int(rows[1].exposure_id)

            custom_id, _ = ensure_custom_group!(db, exposure_id)

            # undoes_event_id: find the most recent index_confirmed for this (group_id, index_id).
            prior = Tables.rowtable(DBInterface.execute(db, """
                SELECT id FROM user_actions
                WHERE action = 'index_confirmed'
                  AND entity_type = 'exposure' AND entity_id = ?
                  AND payload IS NOT NULL
                  AND json_extract(payload, '\$.group_id') = ?
                  AND json_extract(payload, '\$.index_id') = ?
                ORDER BY id DESC LIMIT 1
            """, [exposure_id, custom_id, index_id]))
            undoes = isempty(prior) ? nothing : Int(prior[1].id)

            result = apply_event!(InTransaction(), db, req;
                kind            = "index_unconfirmed",
                entity_type     = "exposure",
                entity_id       = exposure_id,
                payload         = Dict(:group_id => custom_id, :index_id => index_id),
                undoes_event_id = undoes)
            _enqueue_broadcast_from_result!(result, "index_unconfirmed", "exposure", exposure_id)

            body = _group_with_members(db, custom_id)
            body[:event_id]    = result.event_id
            body[:view_row_id] = result.view_row_id
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(body))
        end
    end

    @get "/api/exposures/{id}/speculative-snap" function(req::HTTP.Request, id::Int)
        db = current_db()
        params = HTTP.queryparams(HTTP.URI(req.target))
        phase_name    = get(params, "phase", "")
        anchor_pid_s  = get(params, "anchor_peak_id", "")
        anchor_rp_s   = get(params, "anchor_ratio", "1")

        P = resolve_phase(phase_name)
        P === nothing && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "unknown phase: $phase_name")))

        anchor_peak_id = tryparse(Int, anchor_pid_s)
        anchor_ratio   = something(tryparse(Int, anchor_rp_s), 1)
        anchor_peak_id === nothing && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "missing or invalid anchor_peak_id")))

        peak_rows = Tables.rowtable(DBInterface.execute(db, """
            SELECT a.id, a.q, a.sharpness
            FROM auto_peaks a
            WHERE a.exposure_id = ?
              AND NOT EXISTS (
                  SELECT 1 FROM peak_curations c
                  WHERE c.exposure_id = a.exposure_id AND c.kind = 'exclude'
                    AND ABS(c.q - a.q) <= MAX(1e-6, ABS(a.q) * 0.001)
              )
            UNION ALL
            SELECT id, q, NULL AS sharpness
            FROM peak_curations
            WHERE exposure_id = ? AND kind = 'add'
        """, [id, id]))
        anchor = nothing
        for pr in peak_rows
            if Int(pr.id) == anchor_peak_id
                anchor = pr
                break
            end
        end
        anchor === nothing && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "anchor peak not found in exposure")))

        ratios = Himalaya.phaseratios(P; normalize = true)
        (1 <= anchor_ratio <= length(ratios)) || return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "anchor_ratio out of range for phase")))

        # Exclude the anchor itself from the candidate set so it doesn't
        # snap to a different ratio position than the user requested.
        non_anchor = filter(pr -> Int(pr.id) != anchor_peak_id, peak_rows)
        snaps = compute_snap(non_anchor, P, Float64(anchor.q), anchor_ratio)

        out = map(snaps) do s
            Dict(
                :ratio_position     => s.ratio_position,
                :predicted_q        => s.predicted_q,
                :suggested_peak_id  => s.suggested_peak_id,
                :suggested_q        => s.suggested_q,
                :suggested_residual => s.suggested_residual,
                :is_anchor          => s.ratio_position == anchor_ratio,
            )
        end
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
    end

    @post "/api/exposures/{id}/speculative" function(req::HTTP.Request, id::Int)
        db = current_db()
        return with_idempotency(db, req) do
            body = json(req)

            phase_name      = String(body.phase)
            anchor_peak_id  = Int(body.anchor_peak_id)
            anchor_ratio    = Int(body.anchor_ratio)
            # additional_peak_ids: parallel arrays {ratio_position, peak_id}
            additional      = haskey(body, :additional) ? body.additional : []
            # Default to *not* in the active set — speculative indices are
            # hypotheses, and active membership is an explicit user gesture.
            active_default  = haskey(body, :active) ? Bool(body.active) : false

            P = resolve_phase(phase_name)
            P === nothing && return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "unknown phase: $phase_name")))

            ratio_to_peak = Dict{Int, Int}(anchor_ratio => anchor_peak_id)
            for entry in additional
                rp  = Int(entry.ratio_position)
                pid = Int(entry.peak_id)
                haskey(ratio_to_peak, rp) && return HTTP.Response(400,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "duplicate ratio_position $rp")))
                ratio_to_peak[rp] = pid
            end

            # Run the durable writes inside with_idempotency's outer tx via
            # InTransaction. Defer broadcast on the apply_event! so subscribers
            # don't see uncommitted state if the tx rolls back; queue a
            # post-commit broadcast that fires only after the outer tx commits.
            local nid::Int
            try
                nid = insert_speculative_index!(db, id, P, ratio_to_peak)
            catch e
                return HTTP.Response(400,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => sprint(showerror, e))))
            end
            if active_default
                cid, _ = ensure_custom_group!(db, id)
                DBInterface.execute(db,
                    "INSERT OR IGNORE INTO index_group_members (group_id, index_id) VALUES (?, ?)",
                    [cid, nid])
            end
            payload = Dict(:index_id => nid)
            result = apply_event!(InTransaction(), db, req;
                kind        = "speculative_created",
                entity_type = "exposure",
                entity_id   = id,
                payload     = payload)

            # Defer the SSE broadcast until the outer with_idempotency tx
            # commits. Subscribers converge via cache invalidation in
            # applyRemoteToCache; no post_state needed.
            _enqueue_post_commit_broadcast!(
                Int(result.event_id), "speculative_created", "exposure", Int(id),
                result.user_id, result.client_id, result.client_op_id,
                result.payload_json)

            # Return the freshly-built index in the same shape as GET /api/indices/:id
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT * FROM indices WHERE id = ?", [nid]))
            ix = rows[1]
            peak_rows = Tables.rowtable(DBInterface.execute(db,
                """SELECT ip.peak_id, ip.ratio_position, ip.residual,
                          COALESCE(ap.q, pc.q) AS q_observed
                   FROM index_peaks ip
                   LEFT JOIN auto_peaks ap     ON ap.id = ip.peak_id AND ip.peak_kind = 'auto'
                   LEFT JOIN peak_curations pc ON pc.id = ip.peak_id AND ip.peak_kind = 'curation'
                   WHERE ip.index_id = ? ORDER BY ip.ratio_position""", [nid]))
            predicted = predicted_q_for_phase(String(ix.phase), Float64(ix.basis))
            d = row_to_json(ix)
            d[:peaks]       = rows_to_json(peak_rows)
            d[:predicted_q] = predicted
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(d))
        end
    end

    @delete "/api/indices/{id}" function(req::HTTP.Request, id::Int)
        db = current_db()
        return with_idempotency(db, req) do
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, kind, exposure_id FROM indices WHERE id = ?", [id]))
            isempty(rows) && return HTTP.Response(404,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "index not found")))
            kind_val = String(rows[1].kind)
            kind_val == "speculative" || return HTTP.Response(403,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "only speculative indices can be deleted; use group exclusion for auto indices")))

            # Capture exposure_id BEFORE the DELETE — after deletion the row is gone.
            exposure_id = Int(rows[1].exposure_id)

            DBInterface.execute(db,
                "DELETE FROM index_group_members WHERE index_id = ?", [id])
            DBInterface.execute(db,
                "DELETE FROM index_peaks WHERE index_id = ?", [id])
            DBInterface.execute(db,
                "DELETE FROM indices WHERE id = ?", [id])

            payload = Dict(:index_id => id)
            result = apply_event!(InTransaction(), db, req;
                kind        = "speculative_deleted",
                entity_type = "exposure",
                entity_id   = exposure_id,
                payload     = payload)

            _enqueue_post_commit_broadcast!(
                Int(result.event_id), "speculative_deleted", "exposure",
                Int(exposure_id),
                result.user_id, result.client_id, result.client_op_id,
                result.payload_json)

            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:deleted => id)))
        end
    end

    @get "/api/indices/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM indices WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "index not found")))
        ix        = rows[1]
        peak_rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT ip.peak_id, ip.ratio_position, ip.residual,
                      COALESCE(ap.q, pc.q) AS q_observed
               FROM index_peaks ip
               LEFT JOIN auto_peaks ap     ON ap.id = ip.peak_id AND ip.peak_kind = 'auto'
               LEFT JOIN peak_curations pc ON pc.id = ip.peak_id AND ip.peak_kind = 'curation'
               WHERE ip.index_id = ? ORDER BY ip.ratio_position""", [id]))
        predicted = predicted_q_for_phase(String(ix.phase), Float64(ix.basis))
        d               = row_to_json(ix)
        d[:peaks]       = rows_to_json(peak_rows)
        d[:predicted_q] = predicted
        d[:ngc]         = _ngc_for_phase(String(ix.phase), ix.lattice_d)
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(d))
    end
end
