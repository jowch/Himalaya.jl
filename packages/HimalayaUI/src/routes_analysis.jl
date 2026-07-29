using HTTP, JSON3, DBInterface, Tables, Oxygen, SQLite
using Himalaya

# D-10 (plotting redesign): the legacy index-group machinery — ensure_custom_group!,
# _group_with_members, and the /groups routes — was retired. confirmed_index now
# sources from the durable per-exposure assignment (assignments / assignment_members),
# and the cart uses the assignment-native routes below. The index_groups /
# index_group_members tables are kept (still written by persist_analysis!'s auto
# group, read by migrate_assignments!) but are no longer served or curated.

"""
    _assignment_body(db, exposure_id) -> Dict

Build the canonical assignment response: the durable 3-state assignment for an
exposure plus its 0..N member index ids (ascending). Returns the neutral
default (state 'indexed', no members) when no assignment row exists yet.
"""
function _assignment_body(db::SQLite.DB, exposure_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT state FROM assignments WHERE exposure_id = ?", [exposure_id]))
    state = isempty(rows) ? "indexed" : String(rows[1].state)
    members = Tables.rowtable(DBInterface.execute(db,
        "SELECT index_id FROM assignment_members WHERE exposure_id = ? ORDER BY index_id",
        [exposure_id]))
    Dict(:exposure_id => exposure_id,
         :state       => state,
         :members     => [Int(m.index_id) for m in members])
end

"""
    _assignment_post_state(b) -> Dict

The `{state, members}` envelope carried under `post_state[:assignment]` by the
assignment_* SSE frames. Shared with `_enrich_curation_post_state`
(routes_peaks.jl) so peak-route frames serialize the assignment in exactly the
same shape — one serializer, not two.
"""
_assignment_post_state(b::AbstractDict) =
    Dict(:state => b[:state], :members => b[:members])

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

"""
    _bonnet_for_index(db, exposure_id, phase_name, lattice_d, index_id) -> Union{Dict,Nothing}

For a candidate index, return `{predicted_a, consistent}` when the exposure's
assignment contains a bicontinuous cubic of a DIFFERENT phase that predicts this
candidate's lattice via the Gauss–Bonnet ratio (`Himalaya.bonnet_lattice`).
Returns `nothing` when there is no applicable anchor (no assigned cubic, same
phase, or non-bicontinuous phase). The anchor index itself is never flagged.

Display-and-ranking affordance only — recomputed per request from the live
assignment, never persisted, never folded into `score()` (see docs/scoring.md).
"""
function _bonnet_for_index(db::SQLite.DB, exposure_id::Integer,
                           phase_name::AbstractString, lattice_d, index_id::Integer)
    (lattice_d === nothing || ismissing(lattice_d)) && return nothing
    P = resolve_phase(phase_name)
    P === nothing && return nothing
    a = Float64(lattice_d)
    a > 0 || return nothing

    # Assigned bicontinuous-cubic anchors (phase + lattice), excluding this index.
    anchors = Tables.rowtable(DBInterface.execute(db,
        """SELECT i.phase, i.lattice_d
           FROM assignment_members m JOIN indices i ON i.id = m.index_id
           WHERE m.exposure_id = ? AND i.id != ?
             AND i.lattice_d IS NOT NULL
             AND i.phase IN ('Pn3m', 'Im3m', 'Ia3d')""",
        [Int(exposure_id), Int(index_id)]))

    for anc in anchors
        Pa = resolve_phase(String(anc.phase))
        Pa === nothing && continue
        Pa === P && continue   # same phase: not a coexisting pair
        pred = Himalaya.bonnet_lattice(Pa, Float64(anc.lattice_d), P)
        pred === nothing && continue
        return Dict(:predicted_a => pred,
                    :consistent  => abs(a - pred) <= 0.02 * pred)
    end
    nothing
end

function register_analysis_routes!()
    @get "/api/exposures/{id}/indices" function(req::HTTP.Request, id::Int)
        db = current_db()
        indices = Tables.rowtable(DBInterface.execute(db,
            """SELECT id, exposure_id, phase, basis, score, r_squared, lattice_d,
                  status, kind, inputs_hash
           FROM indices WHERE exposure_id = ? ORDER BY score DESC""", [id]))
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
            d[:bonnet]       = _bonnet_for_index(db, id, String(ix.phase), ix.lattice_d, Int(ix.id))
            d
        end
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(out))
    end

    @get "/api/exposures/{id}/assignment" function(req::HTTP.Request, id::Int)
        db = current_db()
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(_assignment_body(db, id)))
    end

    @post "/api/exposures/{id}/assignment/state" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if !haskey(body, :state)
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "missing field: state")))
        end
        state = String(body.state)
        if !(state in ("indexed", "form_factor", "null"))
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "invalid state: $state")))
        end
        return with_idempotency(db, req) do
            result = apply_event!(InTransaction(), db, req;
                kind        = "assignment_set_state",
                entity_type = "exposure",
                entity_id   = id,
                payload     = Dict(:state => state))
            # Build the response from the now-current assignment, and carry it as
            # the SSE post_state so Plan D's applyRemoteToCache can patch the
            # assignment cache directly (no extra refetch). NOTE: the post_state
            # has NO top-level `indices` key — that is what lets the frontend's
            # CurationPostState cast bail harmlessly for assignment frames.
            b = _assignment_body(db, id)
            post_state = Dict(:assignment => _assignment_post_state(b))
            _enqueue_broadcast_from_result!(result, "assignment_set_state", "exposure", id;
                post_state = post_state)
            b[:event_id]    = result.event_id
            b[:view_row_id] = result.view_row_id
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(b))
        end
    end

    # ── Plan D Task D-1: assignment-native member routes ────────────────────
    # POST/DELETE /api/exposures/{id}/assignment/members are the assignment-
    # native targets that the Plan D frontend calls directly (the legacy
    # /groups/{id}/members dual-write stays live until D-10). They mirror the
    # /assignment/state route shape: with_idempotency + apply_event!(InTransaction())
    # + _assignment_body response, and carry the DISTINCT post_state
    # {assignment:{state,members}} (NO top-level `indices` key) so the
    # frontend's applyPostStateOnly/CurationPostState guard bails harmlessly.
    @post "/api/exposures/{id}/assignment/members" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        if !haskey(body, :index_id)
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "missing field: index_id")))
        end
        local index_id::Int
        try
            index_id = Int(body.index_id)
        catch
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "index_id must be an integer")))
        end
        return with_idempotency(db, req) do
            result = apply_event!(InTransaction(), db, req;
                kind        = "assignment_add",
                entity_type = "exposure",
                entity_id   = id,
                payload     = Dict(:index_id => index_id))
            b = _assignment_body(db, id)
            _enqueue_broadcast_from_result!(result, "assignment_add", "exposure", id;
                post_state = Dict(:assignment => _assignment_post_state(b)))
            b[:event_id]    = result.event_id
            b[:view_row_id] = result.view_row_id
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(b))
        end
    end

    @delete "/api/exposures/{id}/assignment/members/{index_id}" function(req::HTTP.Request, id::Int, index_id::Int)
        db = current_db()
        return with_idempotency(db, req) do
            result = apply_event!(InTransaction(), db, req;
                kind        = "assignment_remove",
                entity_type = "exposure",
                entity_id   = id,
                payload     = Dict(:index_id => index_id))
            b = _assignment_body(db, id)
            _enqueue_broadcast_from_result!(result, "assignment_remove", "exposure", id;
                post_state = Dict(:assignment => _assignment_post_state(b)))
            b[:event_id]    = result.event_id
            b[:view_row_id] = result.view_row_id
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(b))
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
        body = json(req)
        for field in (:phase, :anchor_peak_id, :anchor_ratio)
            if !haskey(body, field)
                return HTTP.Response(400,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "missing field: $(field)")))
            end
        end
        local phase_name::String
        local anchor_peak_id::Int
        local anchor_ratio::Int
        try
            phase_name     = String(body.phase)
            anchor_peak_id = Int(body.anchor_peak_id)
            anchor_ratio   = Int(body.anchor_ratio)
        catch
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "phase must be string; anchor_peak_id and anchor_ratio must be integers")))
        end
        return with_idempotency(db, req) do
            # additional_peak_ids: parallel arrays {ratio_position, peak_id}
            additional      = haskey(body, :additional) ? body.additional : []
            # NOTE: the `active` body field is accepted for back-compat but no
            # longer auto-adds to any set — D-10 retired the legacy active custom
            # group, and "make this index active" is now an explicit assignment
            # gesture (POST /assignment/members). The only caller already passes
            # active:false, so this is a no-op contract change.

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
                """SELECT id, exposure_id, phase, basis, score, r_squared, lattice_d,
                  status, kind, inputs_hash
           FROM indices WHERE id = ?""", [nid]))
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
            # Match the shape of GET /api/indices/:id (which sets ngc on the
            # response). Without this, `createSpeculativeMutator.onSuccess`
            # writes an IndexEntry-without-`ngc` into the indices cache, and
            # downstream `index.ngc` reads return undefined for any
            # speculatively-created index until a refetch.
            d[:ngc]         = _ngc_for_phase(String(ix.phase), ix.lattice_d)
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(d))
        end
    end

    # ── Plan D Task D-9 (B4): client-fitted custom-index commit ─────────────
    # Accepts {phase, basis, ratios?} where `basis` is the q₁ slope the modal
    # computed from a symmetry + lattice via real physics, and `ratios` are the
    # NORMALIZED ratios it DREW. Persists a basis-locked speculative index,
    # claims the peaks landing within CUSTOM_SNAP_TOL of those reflections, and
    # adds it to the assignment in one transaction. A ratio set rather than a
    # count because the modal's series is not a positional prefix of the core
    # one for Hexagonal (see compute_snap). Two SSE frames: speculative_created (indices
    # cache) THEN assignment_add (assignment cache). ORDERING IS LOAD-BEARING —
    # the own-tab deferred resolves off the FIRST frame (speculative_created)
    # which carries the new IndexEntry, then assignment_add patches the cart.
    @post "/api/exposures/{id}/custom-index" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        for field in (:phase, :basis)
            if !haskey(body, field)
                return HTTP.Response(400, ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "missing field: $(field)")))
            end
        end
        local phase_name::String
        local basis::Float64
        local ratios::Union{Vector{Float64}, Nothing}
        try
            phase_name = String(body.phase)
            basis      = Float64(body.basis)
            # Optional for back-compat: absent means "scan the whole ratio
            # series", which is only correct for a caller with no truncated
            # display of its own. The modal always sends it.
            # `isa AbstractVector` is load-bearing, not defensive: a JSON
            # STRING iterates as Chars and Float64('s') == 115.0, so "six"
            # would convert to three positive "ratios", clear every check
            # below, match zero positions in compute_snap, and commit an index
            # claiming NO peaks — a silent 200 reintroducing the exact bug this
            # route exists to fix.
            raw = haskey(body, :ratios) ? body.ratios : nothing
            raw === nothing || raw isa AbstractVector ||
                error("ratios must be an array")
            ratios     = raw === nothing ? nothing :
                         Float64[Float64(r) for r in raw]
        catch
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "phase must be string; basis a number; ratios an array of numbers")))
        end
        if ratios !== nothing && (isempty(ratios) || any(r -> !(r > 0), ratios))
            return HTTP.Response(400, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "ratios must be a non-empty array of positive numbers")))
        end
        P = resolve_phase(phase_name)
        P === nothing && return HTTP.Response(400, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "unknown phase: $phase_name")))
        basis > 0 || return HTTP.Response(400, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "basis must be positive")))

        return with_idempotency(db, req) do
            # No catch → 400 here: this body runs INSIDE with_idempotency's
            # transaction, and a normal return COMMITS (idempotency.jl only
            # skips the cache row on >=400, it does not roll back). Now that
            # insert_custom_index! issues 1+2N statements, swallowing a
            # mid-loop throw would commit an orphan `indices` row with no
            # event. Let it propagate so the transaction unwinds; the
            # validation above already rejects the reachable bad inputs.
            nid = insert_custom_index!(db, id, P, basis; drawn_ratios = ratios)

            # Frame 1: speculative_created — carries the new index to the indices
            # cache (post-commit broadcast; subscribers converge via invalidate).
            sc = apply_event!(InTransaction(), db, req;
                kind        = "speculative_created",
                entity_type = "exposure",
                entity_id   = id,
                payload     = Dict(:index_id => nid))
            _enqueue_post_commit_broadcast!(
                Int(sc.event_id), "speculative_created", "exposure", Int(id),
                sc.user_id, sc.client_id, sc.client_op_id, sc.payload_json)

            # Frame 2: assignment_add — adds the new custom index to the cart,
            # carrying the distinct {assignment} post_state.
            aa = apply_event!(InTransaction(), db, req;
                kind        = "assignment_add",
                entity_type = "exposure",
                entity_id   = id,
                payload     = Dict(:index_id => nid))
            ab = _assignment_body(db, id)
            _enqueue_broadcast_from_result!(aa, "assignment_add", "exposure", id;
                post_state = Dict(:assignment =>
                    Dict(:state => ab[:state], :members => ab[:members])))

            # Return the freshly-built index (GET /api/indices/:id shape).
            rows = Tables.rowtable(DBInterface.execute(db,
                """SELECT id, exposure_id, phase, basis, score, r_squared, lattice_d,
                      status, kind, inputs_hash
               FROM indices WHERE id = ?""", [nid]))
            ix = rows[1]
            # insert_custom_index! claims the peaks the modal showed landing —
            # read them back rather than asserting []; a hardcoded empty list
            # left the Focus comb/detector/cart showing zero matches.
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
            d[:ngc]         = _ngc_for_phase(String(ix.phase), ix.lattice_d)
            # Match GET /api/exposures/:id/indices, which sets `bonnet` (:157).
            # CandidateRow reads `ix.bonnet?.consistent`, so omitting it here
            # left a freshly committed custom index without its Gauss–Bonnet
            # flag until the next refetch — the same hazard this file already
            # documents for `ngc` on the sibling speculative route.
            d[:bonnet]      = _bonnet_for_index(db, id, String(ix.phase), ix.lattice_d, nid)
            d[:event_id]    = sc.event_id
            d[:view_row_id] = sc.view_row_id
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
            # Likewise capture membership before the cascade wipes it: the
            # assignment_remove frame below is only honest if the index really
            # was in the call.
            was_member = !isempty(Tables.rowtable(DBInterface.execute(db,
                "SELECT 1 FROM assignment_members WHERE exposure_id = ? AND index_id = ?",
                [exposure_id, id])))

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

            # assignment_members.index_id is ON DELETE CASCADE, so the DELETE
            # above silently dropped the membership row. events.jl declares the
            # assignment_* kinds the SOLE writer to that table and
            # speculative_deleted is a dispatcher no-op, so without this frame
            # rebuild_views_from_log! would re-fold a member row whose `indices`
            # row no longer exists (INSERT OR IGNORE does not suppress an FK
            # violation). Emitted AFTER speculative_deleted so the own-tab
            # deferred still resolves off the first frame, matching the
            # custom-index route's ordering contract, and it carries the
            # canonical {state, members} post_state so subscribers patch rather
            # than invalidate.
            if was_member
                ar = apply_event!(InTransaction(), db, req;
                    kind        = "assignment_remove",
                    entity_type = "exposure",
                    entity_id   = exposure_id,
                    payload     = Dict(:index_id => id))
                ab = _assignment_body(db, exposure_id)
                _enqueue_broadcast_from_result!(ar, "assignment_remove", "exposure", exposure_id;
                    post_state = Dict(:assignment => _assignment_post_state(ab)))
            end

            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:deleted => id)))
        end
    end

    @get "/api/indices/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT id, exposure_id, phase, basis, score, r_squared, lattice_d,
                  status, kind, inputs_hash
           FROM indices WHERE id = ?""", [id]))
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
