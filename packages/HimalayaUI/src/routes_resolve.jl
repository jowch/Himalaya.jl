using HTTP, JSON3, DBInterface, Tables, Oxygen

# ─────────────────────────────────────────────────────────────────────────────
# Spec §3.1 — `/api/resolve` slug-to-id resolver.
#
# Read-only. No writes; no with_idempotency; no apply_event!; no SSE; no
# client_op_id; no pendingDeferreds participation. Three SELECTs and a
# response. Queue-orthogonal by design — confirmed against
# docs/mutation-queue.md. Future maintainers tempted to extend it with a
# write path: please don't; add a sibling endpoint instead.
#
# Tiebreaker for non-unique experiment names: ORDER BY id ASC LIMIT 1
# (deterministic across runs). samples and exposures are uniqueness-
# enforced post-#88 / Task 1.
# ─────────────────────────────────────────────────────────────────────────────

function _json(status::Int, body)
    HTTP.Response(status, ["Content-Type" => "application/json"], JSON3.write(body))
end

function _has_param(params, name::String)
    haskey(params, name) && !isempty(params[name])
end

function _safe_str(v)
    # Tables.rowtable returns `missing` for SQL NULL. Coerce to "" so the
    # response shape stays string-typed; the caller may post-process.
    ismissing(v) ? "" : String(v)
end

function _parse_id_or_400(s::String, field::String)
    n = tryparse(Int, s)
    n === nothing && return _json(400, Dict(:error => "invalid_id", :field => field))
    return n
end

function register_resolve_routes!()
    @get "/api/resolve" function(req::HTTP.Request)
        db = current_db()
        params = HTTP.queryparams(HTTP.URI(req.target))

        # Mutual-exclusion check: same-entity name+id collision is 400.
        for (n, i) in (("experiment", "experiment_id"),
                       ("sample",     "sample_id"),
                       ("exposure",   "exposure_id"))
            if _has_param(params, n) && _has_param(params, i)
                return _json(400, Dict(:error => "ambiguous_params"))
            end
        end

        # Resolve experiment. NULL-name rows are treated as "no canonical
        # slug" → 404. The frontend has no path to construct a URL for an
        # experiment without a name; rejecting them here keeps the round
        # trip well-formed.
        exp_row = nothing
        if _has_param(params, "experiment")
            name = params["experiment"]
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, name FROM experiments WHERE name = ? ORDER BY id ASC LIMIT 1", [name]))
            isempty(rows) && return _json(404, Dict(
                :error => "not_found", :missing => "experiment", :missing_value => name))
            exp_row = (id=Int(rows[1].id), name=_safe_str(rows[1].name))
        elseif _has_param(params, "experiment_id")
            id_or_resp = _parse_id_or_400(params["experiment_id"], "experiment_id")
            id_or_resp isa HTTP.Response && return id_or_resp
            id = id_or_resp::Int
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, name FROM experiments WHERE id = ? LIMIT 1", [id]))
            isempty(rows) && return _json(404, Dict(
                :error => "not_found", :missing => "experiment",
                :missing_value => string(id)))
            ismissing(rows[1].name) && return _json(404, Dict(
                :error => "not_found", :missing => "experiment",
                :missing_value => string(id),
                :reason => "experiment has no canonical name"))
            exp_row = (id=Int(rows[1].id), name=_safe_str(rows[1].name))
        else
            return _json(400, Dict(:error => "missing_experiment"))
        end

        # Resolve sample (optional).
        sample_row = nothing
        if _has_param(params, "sample") || _has_param(params, "sample_id")
            if _has_param(params, "sample")
                name = params["sample"]
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, name FROM samples WHERE experiment_id = ? AND name = ? LIMIT 1",
                    [exp_row.id, name]))
                isempty(rows) && return _json(404, Dict(
                    :error => "not_found", :missing => "sample", :missing_value => name,
                    :experiment_resolved => Dict(:id => exp_row.id, :name => exp_row.name)))
                sample_row = (id=Int(rows[1].id), name=_safe_str(rows[1].name))
            else
                id_or_resp = _parse_id_or_400(params["sample_id"], "sample_id")
                id_or_resp isa HTTP.Response && return id_or_resp
                id = id_or_resp::Int
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, name FROM samples WHERE id = ? AND experiment_id = ? LIMIT 1",
                    [id, exp_row.id]))
                isempty(rows) && return _json(404, Dict(
                    :error => "not_found", :missing => "sample", :missing_value => string(id),
                    :experiment_resolved => Dict(:id => exp_row.id, :name => exp_row.name)))
                sample_row = (id=Int(rows[1].id), name=_safe_str(rows[1].name))
            end
        end

        # Resolve exposure (optional, requires sample).
        exposure_row = nothing
        if _has_param(params, "exposure") || _has_param(params, "exposure_id")
            sample_row === nothing && return _json(400, Dict(:error => "missing_sample"))
            if _has_param(params, "exposure")
                name = params["exposure"]
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, filename FROM exposures WHERE sample_id = ? AND filename = ? LIMIT 1",
                    [sample_row.id, name]))
                isempty(rows) && return _json(404, Dict(
                    :error => "not_found", :missing => "exposure", :missing_value => name,
                    :experiment_resolved => Dict(:id => exp_row.id, :name => exp_row.name),
                    :sample_resolved     => Dict(:id => sample_row.id, :name => sample_row.name)))
                exposure_row = (id=Int(rows[1].id), filename=_safe_str(rows[1].filename))
            else
                id_or_resp = _parse_id_or_400(params["exposure_id"], "exposure_id")
                id_or_resp isa HTTP.Response && return id_or_resp
                id = id_or_resp::Int
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, filename FROM exposures WHERE id = ? AND sample_id = ? LIMIT 1",
                    [id, sample_row.id]))
                isempty(rows) && return _json(404, Dict(
                    :error => "not_found", :missing => "exposure", :missing_value => string(id),
                    :experiment_resolved => Dict(:id => exp_row.id, :name => exp_row.name),
                    :sample_resolved     => Dict(:id => sample_row.id, :name => sample_row.name)))
                exposure_row = (id=Int(rows[1].id), filename=_safe_str(rows[1].filename))
            end
        end

        # Build success response.
        body = Dict{Symbol,Any}(
            :experiment_id   => exp_row.id,
            :experiment_name => exp_row.name,
        )
        if sample_row !== nothing
            body[:sample_id]   = sample_row.id
            body[:sample_name] = sample_row.name
        end
        if exposure_row !== nothing
            body[:exposure_id]       = exposure_row.id
            body[:exposure_filename] = exposure_row.filename
        end
        _json(200, body)
    end
end
