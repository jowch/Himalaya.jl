using HTTP, JSON3, DBInterface, Tables, Oxygen, SQLite

function _resolve_analysis_dir(db::SQLite.DB, exposure_id::Int)::Union{String, Nothing}
    setup = exposure_setup(db, exposure_id)
    setup === nothing ? nothing : String(setup.analysis_dir)
end

"""
    _enrich_curation_post_state(db, exposure_id;
                                dropped_assignment_phases = String[]) -> Dict

Read the post-mutation state needed by the M2 SSE post_state contract:
the current `analysis_inputs_hash` for the exposure plus the current
`indices` array, and (F-WIPE W1) the current `assignment` — in the same
`{state, members}` shape the assignment_* frames carry (shared
`_assignment_post_state` serializer) — so foreign tabs converge on the
re-attached membership without a refetch. `assignment_dropped` is included
ONLY when the reanalysis dropped members (phases whose semantic identity
failed to re-attach; threaded from `analyze_exposure!`'s return). All reads
happen inside the caller's transaction so the snapshot is consistent.
"""
function _enrich_curation_post_state(db::SQLite.DB, exposure_id::Int;
                                     dropped_assignment_phases::Vector{String} = String[])
    ps = Dict{Symbol, Any}(
        :analysis_inputs_hash => read_inputs_hash(db, exposure_id),
        :indices              => _serialized_indices_for_broadcast(db, exposure_id),
        :assignment           => _assignment_post_state(_assignment_body(db, exposure_id)),
    )
    isempty(dropped_assignment_phases) ||
        (ps[:assignment_dropped] = dropped_assignment_phases)
    ps
end

"""
    _queue_peak_added_broadcast!(result, entity_id, payload, post_state)

Specialized broadcast helper for `peak_added`. The payload is mutated in the
route body to include `peak_curation_id` AFTER `apply_event!` returns (the
id is the dispatcher's `view_row_id`), so we have to re-serialize the dict
here — `result.payload_json` was frozen by `apply_event!` before the
mutation. All other curation kinds use `_enqueue_broadcast_from_result!`
directly (it reuses `result.payload_json` with no re-serialization —
suggestion #4 from PR review).
"""
function _queue_peak_added_broadcast!(result, entity_id::Integer,
                                       payload::AbstractDict, post_state)
    _enqueue_post_commit_broadcast!(
        Int(result.event_id), "peak_added", "exposure", Int(entity_id),
        result.user_id, result.client_id, result.client_op_id,
        JSON3.write(payload);
        post_state = post_state)
end

function _auto_peak_row(db::SQLite.DB, id::Int)
    owner = Tables.rowtable(DBInterface.execute(db,
        "SELECT exposure_id FROM auto_peaks WHERE id = ?", [id]))
    isempty(owner) && return nothing
    rows = peak_listing(db, Int(owner[1].exposure_id); include_superseded = true)
    i = findfirst(r -> r.source == "auto" && Int(r.id) == id, rows)
    i === nothing ? nothing : rows[i]
end

function register_peaks_routes!()
    @get "/api/exposures/{id}/peaks" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = peak_listing(db, id)
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(rows_to_json(rows; bool_keys = (:excluded,))))
    end

    @post "/api/exposures/{id}/peaks" function(req::HTTP.Request, id::Int)
        db = current_db()
        body = json(req)
        # Validate before with_idempotency so a malformed body returns 400
        # (validation toast on the frontend) rather than 500 (infrastructure
        # banner). The frontend's failure-class router treats 4xx as
        # validation; an uncaught Float64() conversion throw would otherwise
        # surface as 500 and mis-route. Caught by the smoke checklist.
        if !haskey(body, :q)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "missing field: q")))
        end
        local q::Float64
        try
            q = Float64(body.q)
        catch
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "q must be a number")))
        end
        return with_idempotency(db, req) do

            # Initial payload before insertion. The peak_curation_id is added
            # below once the dispatcher returns the new row id, so the SSE
            # frame is self-contained — foreign tabs can use the real id
            # rather than synthesizing a placeholder from event_id (issue #2).
            payload = Dict{Symbol, Any}(:q => q)

            # Participate in with_idempotency's outer transaction so the event
            # write, the curation insert, the synchronous reanalyze, AND the
            # idempotent_responses cache row all commit (or roll back) together.
            result = apply_event!(InTransaction(), db, req;
                kind        = "peak_added",
                entity_type = "exposure",
                entity_id   = id,
                payload     = payload)

            new_curation_id = result.view_row_id
            new_curation_id === nothing &&
                error("dispatcher did not record a view_row_id for peak_added")
            payload[:peak_curation_id] = new_curation_id

            # Synchronous reanalyze inside the tx; trace is unchanged so use
            # the fast-skip path. defer_broadcast=true suppresses the inner
            # analyze_run frame — the curation broadcast carries the post_state.
            dir = _resolve_analysis_dir(db, id)
            dropped = String[]
            if dir !== nothing
                ares = analyze_exposure!(db, id, dir;
                    trace_known_unchanged = true,
                    defer_broadcast       = true)
                dropped = ares.dropped_assignment_phases
            end

            new_hash  = read_inputs_hash(db, id)
            post_state = _enrich_curation_post_state(db, id;
                dropped_assignment_phases = dropped)

            _queue_peak_added_broadcast!(result, id, payload, post_state)

            HTTP.Response(201, ["Content-Type" => "application/json"],
                JSON3.write(Dict(
                    # Full Peak-shape — manual peaks have no analysis-derived
                    # metrics yet, so intensity/prominence/sharpness are
                    # explicitly null (NOT omitted). Pre-fix the response
                    # omitted these three; PeakAddResponse extends Peak so
                    # the cache write would leave them as `undefined`,
                    # silently violating the type's `null` contract.
                    # Caught in sixth-pass review (issue #17).
                    :id                   => new_curation_id,
                    :exposure_id          => id,
                    :q                    => q,
                    :intensity            => nothing,
                    :prominence           => nothing,
                    :sharpness            => nothing,
                    :source               => "manual",
                    :excluded             => false,
                    :event_id             => result.event_id,
                    :view_row_id          => new_curation_id,
                    :analysis_inputs_hash => new_hash,
                )))
        end
    end

    @patch "/api/peaks/auto/{id}" function(req::HTTP.Request, id::Int)
        db = current_db()
        return with_idempotency(db, req) do
            body = json(req)

            if !haskey(body, :excluded)
                return HTTP.Response(400,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "missing field: excluded")))
            end

            peak_row = _auto_peak_row(db, id)
            if peak_row === nothing
                return HTTP.Response(404,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "auto peak not found")))
            end

            exposure_id = Int(peak_row.exposure_id)

            excluded = Bool(body.excluded)
            peak_q   = Float64(peak_row.q)
            q_step   = exposure_q_step(db, exposure_id)

            event_id_out::Union{Int, Nothing} = nothing
            view_row_id_out::Union{Int, Nothing} = nothing
            mutated = false

            if excluded
                # Idempotent: only emit an event if a curation already hides this peak.
                existing = Int(peak_row.excluded) == 1 ? [1] : Int[]
                if isempty(existing)
                    payload = Dict(:q => peak_q, :auto_peak_id => id)
                    result = apply_event!(InTransaction(), db, req;
                        kind        = "peak_excluded",
                        entity_type = "exposure",
                        entity_id   = exposure_id,
                        payload     = payload)
                    event_id_out    = result.event_id
                    view_row_id_out = result.view_row_id
                    mutated = true

                    dir = _resolve_analysis_dir(db, exposure_id)
                    dropped = String[]
                    if dir !== nothing
                        ares = analyze_exposure!(db, exposure_id, dir;
                            trace_known_unchanged = true,
                            defer_broadcast       = true)
                        dropped = ares.dropped_assignment_phases
                    end

                    post_state = _enrich_curation_post_state(db, exposure_id;
                        dropped_assignment_phases = dropped)
                    _enqueue_broadcast_from_result!(result, "peak_excluded",
                        "exposure", exposure_id; post_state = post_state)
                end
            elseif Int(peak_row.excluded) == 1
                # peak_unexcluded: find the prior peak_excluded event to record as undoes.
                prior = Tables.rowtable(DBInterface.execute(db, """
                    SELECT id FROM user_actions
                    WHERE action = 'peak_excluded'
                      AND entity_type = 'exposure' AND entity_id = ?
                      AND payload IS NOT NULL
                      AND ABS(json_extract(payload, '\$.q') - ?) <= ?
                    ORDER BY id DESC LIMIT 1
                """, [exposure_id, peak_q, peak_match_tol(peak_q, q_step)]))
                undoes = isempty(prior) ? nothing : Int(prior[1].id)

                payload = Dict(:q => peak_q, :auto_peak_id => id,
                               :curation_qs => exclude_qs_for_peak(db, exposure_id, peak_q))
                result = apply_event!(InTransaction(), db, req;
                    kind            = "peak_unexcluded",
                    entity_type     = "exposure",
                    entity_id       = exposure_id,
                    payload         = payload,
                    undoes_event_id = undoes)
                event_id_out    = result.event_id
                view_row_id_out = result.view_row_id
                mutated = true

                dir = _resolve_analysis_dir(db, exposure_id)
                dropped = String[]
                if dir !== nothing
                    ares = analyze_exposure!(db, exposure_id, dir;
                        trace_known_unchanged = true,
                        defer_broadcast       = true)
                    dropped = ares.dropped_assignment_phases
                end

                post_state = _enrich_curation_post_state(db, exposure_id;
                    dropped_assignment_phases = dropped)
                _enqueue_broadcast_from_result!(result, "peak_unexcluded",
                    "exposure", exposure_id; post_state = post_state)
            end

            # Return the updated row with fresh excluded state, plus the M2
            # response shape (event_id, view_row_id, analysis_inputs_hash).
            out = _auto_peak_row(db, id)
            out === nothing && return HTTP.Response(404,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "auto peak not found after update")))

            d = row_to_json(out; bool_keys = (:excluded,))
            d[:event_id]             = event_id_out
            d[:view_row_id]          = view_row_id_out
            d[:analysis_inputs_hash] = read_inputs_hash(db, exposure_id)
            HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(d))
        end
    end

    @delete "/api/peaks/manual/{id}" function(req::HTTP.Request, id::Int)
        db = current_db()
        return with_idempotency(db, req) do
            curation_rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT exposure_id, q FROM peak_curations WHERE id = ? AND kind = 'add'",
                [id]))
            if isempty(curation_rows)
                return HTTP.Response(404,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "manual peak not found")))
            end
            exposure_id = Int(curation_rows[1].exposure_id)
            removed_q   = Float64(curation_rows[1].q)

            DBInterface.execute(db,
                "DELETE FROM index_peaks WHERE peak_id = ? AND peak_kind = 'curation'", [id])
            DBInterface.execute(db,
                "DELETE FROM peak_curations WHERE id = ?", [id])

            payload = Dict{Symbol, Any}(:peak_curation_id => id, :q => removed_q)
            result = apply_event!(InTransaction(), db, req;
                kind        = "peak_removed",
                entity_type = "exposure",
                entity_id   = exposure_id,
                payload     = payload)

            dir = _resolve_analysis_dir(db, exposure_id)
            dropped = String[]
            if dir !== nothing
                ares = analyze_exposure!(db, exposure_id, dir;
                    trace_known_unchanged = true,
                    defer_broadcast       = true)
                dropped = ares.dropped_assignment_phases
            end

            post_state = _enrich_curation_post_state(db, exposure_id;
                dropped_assignment_phases = dropped)
            _enqueue_broadcast_from_result!(result, "peak_removed",
                "exposure", exposure_id; post_state = post_state)

            HTTP.Response(200,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(
                    :event_id             => result.event_id,
                    :view_row_id          => result.view_row_id,
                    :analysis_inputs_hash => post_state[:analysis_inputs_hash],
                )))
        end
    end
end
