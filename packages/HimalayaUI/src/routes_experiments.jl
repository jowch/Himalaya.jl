using HTTP, JSON3, DBInterface, Tables, Oxygen

const _READONLY_FIELDS = ["data_dir", "analysis_dir", "manifest_path", "path", "id", "created_at"]
const _PLAIN_PATCH_FIELDS = ["name", "description"]
const _CLEARABLE_PATCH_FIELDS = ("description",)

"""
    _experiment_stats(db, exp_id) -> NamedTuple

Cheap roll-up of counts for the shared experiment header stat ledger and the
gallery cards. `span_hours` is the exposure timestamp span; `sessions` is the
persisted macro-session count from `loads.session_id` (assigned at ingest by
`_assign_sessions` in grouping.jl, backfilled for existing loads by
`backfill_load_sessions!` in db.jl — spec §5). `started_at` is the minimum
exposure timestamp (String or nothing) — used by the gallery for year-grouping.
"""
function _experiment_stats(db::SQLite.DB, exp_id::Integer)
    loads = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM loads WHERE experiment_id = ?", [exp_id]))[1].c
    samples = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM samples WHERE experiment_id = ?", [exp_id]))[1].c
    exposures = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM exposures WHERE experiment_id = ?", [exp_id]))[1].c
    span = Tables.rowtable(DBInterface.execute(db,
        "SELECT (julianday(MAX(timestamp)) - julianday(MIN(timestamp))) * 24 AS h " *
        "FROM exposures WHERE experiment_id = ? AND timestamp IS NOT NULL", [exp_id]))[1].h
    span_hours = (span === missing || span === nothing) ? 0.0 : round(Float64(span); digits = 1)
    sessions = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(DISTINCT session_id) AS c FROM loads WHERE experiment_id = ? AND session_id IS NOT NULL",
        [exp_id]))[1].c
    t = Tables.rowtable(DBInterface.execute(db,
        "SELECT MIN(timestamp) AS t FROM exposures WHERE experiment_id = ? AND timestamp IS NOT NULL",
        [exp_id]))[1].t
    started_at = (t === missing || t === nothing) ? nothing : String(t)
    (loads = Int(loads), samples = Int(samples), exposures = Int(exposures),
     sessions = Int(sessions), span_hours = span_hours, started_at = started_at)
end

function _setup_row_to_json(row::NamedTuple)
    d = row_to_json(row)
    d[:included] = Int(row.included) == 1
    d[:name] = setup_display_name(row)
    if ismissing(row.q_units)
        d[:q_units] = "A-1"
        d[:q_units_source] = "default"
    end
    d
end

function _experiment_row_to_json(row::NamedTuple, db::Union{SQLite.DB, Nothing} = nothing)
    d = row_to_json(row)
    if db !== nothing
        exp_id = Int(row.id)
        d[:setups] = [_setup_row_to_json(st) for st in setups_of(db, exp_id)]
        d[:stats] = _experiment_stats(db, exp_id)
    end
    d
end

_body_get(body, k) = get(body, k, get(body, string(k), nothing))

_json_error(status, d) = HTTP.Response(status, ["Content-Type" => "application/json"], JSON3.write(d))

function _setup_spec(spec; default_analysis_dir)
    pats = something(_body_get(spec, :patterns), Dict())
    ppat(k) = (v = _body_get(pats, k); v === nothing ? nothing : String(v))
    geo = something(_body_get(spec, :geometry), Dict())
    gnum(k) = (v = _body_get(geo, k); v isa Number ? Float64(v) : nothing)
    gstr(k) = (v = _body_get(geo, k); v isa AbstractString ? String(v) : nothing)
    psrc(val, k) = val === nothing ? "default" : something(gstr(Symbol(k, :_source)), "user")
    ad = _body_get(spec, :analysis_dir)
    analysis_dir = ad !== nothing ? String(ad) : default_analysis_dir
    included = _body_get(spec, :included)
    out = Dict{Symbol,Any}(:analysis_dir => analysis_dir,
        :image_pattern => ppat(:image), :metadata_pattern => ppat(:metadata),
        :integration_pattern => ppat(:integration),
        :included => included === nothing ? true : included === true)
    for k in (:beam_center_x, :beam_center_y, :flight_path_m, :pixel_size_um, :energy_kev)
        v = gnum(k)
        out[k] = v
        out[Symbol(k, :_source)] = psrc(v, k)
    end
    out
end

_is_obj(v) = v isa AbstractDict || v isa JSON3.Object
_under(child, parent) = (c = _canonical_dir(child); q = _canonical_dir(parent);
    c == q || startswith(c, rstrip(q, '/') * "/"))

function _malformed_setup(s)
    _is_obj(s) || return "each setup must be an object"
    _body_get(s, :data_dir) isa AbstractString || return "each setup needs a data_dir"
    for k in (:root, :analysis_dir)
        v = _body_get(s, k)
        v === nothing || v isa AbstractString || return "$k must be a string"
    end
    inc = _body_get(s, :included)
    inc === nothing || inc isa Bool || return "included must be true or false"
    for k in (:patterns, :geometry)
        v = _body_get(s, k)
        v === nothing || _is_obj(v) || return "$k must be an object"
    end
    _patterns_problem(_body_get(s, :patterns))
end

function _patterns_problem(pats)
    pats === nothing && return nothing
    _is_obj(pats) || return "patterns must be an object"
    for k in (:image, :metadata, :integration)
        bad = pattern_problem(String(k), _body_get(pats, k))
        bad === nothing || return bad
    end
    nothing
end

_default_analysis_dir(data_dir) = (ad = joinpath(data_dir, "analysis"); isdir(ad) ? ad : data_dir)

function _checked_data_dir(db, raw)
    raw === nothing && return _json_error(400, Dict(:error => "data_dir (or path) is required"))
    data_dir = String(raw)
    isdir(data_dir) || return _json_error(400,
        Dict(:error => "path does not exist or is not a directory", :path => data_dir))
    data_dir = _canonical_dir(data_dir)
    existing = experiment_using_dir(db, data_dir)
    existing === nothing || return _json_error(409,
        Dict(:error => "an experiment already uses this directory",
             :path => data_dir, :experiment_id => existing.id))
    data_dir
end

function _create_one_setup_experiment(db, body)
    bad = _patterns_problem(_body_get(body, :patterns))
    bad === nothing || return _json_error(400, Dict(:error => bad))
    data_dir = _checked_data_dir(db, something(_body_get(body, :data_dir), _body_get(body, :path), Some(nothing)))
    data_dir isa HTTP.Response && return data_dir
    name_val = _body_get(body, :name)
    exp_name = name_val !== nothing ? String(name_val) : basename(rstrip(data_dir, '/'))
    spec = _setup_spec(body; default_analysis_dir = _default_analysis_dir(data_dir))
    delete!(spec, :included)
    exp_id = lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            create_experiment!(db; name = exp_name, path = setup_root_dir(data_dir),
                data_dir, ingest_status = "scanning", spec...)
        end
    end
    (exp_id, exp_name, data_dir)
end

function _beamtime_setup_specs(db, path, raw_setups)
    raw_setups isa AbstractVector && !isempty(raw_setups) ||
        return _json_error(400, Dict(:error => "setups must be a non-empty list"))
    specs = Dict{Symbol,Any}[]
    for s in raw_setups
        bad = _malformed_setup(s)
        bad === nothing || return _json_error(400, Dict(:error => bad))
        data_dir = _checked_data_dir(db, _body_get(s, :data_dir))
        data_dir isa HTTP.Response && return data_dir
        any(p -> _same_dir(p[:data_dir], data_dir), specs) && return _json_error(400,
            Dict(:error => "two setups use the same data folder", :path => data_dir))
        spec = _setup_spec(s; default_analysis_dir = _default_analysis_dir(data_dir))
        spec[:data_dir] = data_dir
        root = _body_get(s, :root)
        root_dir = root === nothing ? setup_root_dir(data_dir) : _canonical_dir(String(root))
        isdir(root_dir) && _under(data_dir, root_dir) && _under(root_dir, path) ||
            return _json_error(400, Dict(:error => "a setup folder must hold its data folder and sit in the beamtime folder",
                                         :path => root_dir))
        spec[:root_dir] = root_dir
        push!(specs, spec)
    end
    any(sp -> sp[:included], specs) ||
        return _json_error(400, Dict(:error => "tick at least one setup"))
    specs
end

function _add_beamtime_setups(db, id, body)
    rows = Tables.rowtable(DBInterface.execute(db, "SELECT path FROM experiments WHERE id = ?", [id]))
    isempty(rows) && return _json_error(404, Dict(:error => "experiment not found"))
    path = _canonical_dir(String(first(rows).path))
    any(st -> _same_dir(st.root_dir, path), setups_of(db, id)) && return _json_error(400,
        Dict(:error => "this experiment is one setup folder, so it cannot hold more setups"))
    specs = _beamtime_setup_specs(db, path, _body_get(body, :setups))
    specs isa HTTP.Response && return specs
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            Int[create_setup!(db; experiment_id = id, position = _next_position(db, id), spec...)
                for spec in specs]
        end
    end
end

function _create_beamtime_experiment(db, body, raw_setups)
    raw_path = _body_get(body, :path)
    raw_path === nothing && return _json_error(400, Dict(:error => "path is required"))
    isdir(String(raw_path)) || return _json_error(400,
        Dict(:error => "path does not exist or is not a directory", :path => String(raw_path)))
    path = _canonical_dir(String(raw_path))
    existing = experiment_at_path(db, path)
    existing === nothing || return _json_error(409, Dict(
        :error => "this beamtime is already added as $(existing.name)",
        :experiment_id => existing.id, :name => existing.name, :path => path))
    name_val = _body_get(body, :name)
    name_val === nothing || name_val isa AbstractString ||
        return _json_error(400, Dict(:error => "name must be a string"))
    specs = _beamtime_setup_specs(db, path, raw_setups)
    specs isa HTTP.Response && return specs
    exp_name = name_val !== nothing && !isempty(strip(String(name_val))) ?
        String(strip(String(name_val))) : basename(rstrip(path, '/'))
    exp_id = lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            eid = create_experiment!(db; name = exp_name, path, ingest_status = "scanning")
            for (i, spec) in enumerate(specs)
                create_setup!(db; experiment_id = eid, position = i - 1, spec...)
            end
            eid
        end
    end
    (exp_id, exp_name, specs[1][:data_dir])
end

function register_experiments_routes!()
    @post "/api/experiments" function(req::HTTP.Request)
        db   = current_db()
        body = json(req)
        raw_setups = _body_get(body, :setups)
        created = raw_setups === nothing ? _create_one_setup_experiment(db, body) :
                                           _create_beamtime_experiment(db, body, raw_setups)
        created isa HTTP.Response && return created
        exp_id, exp_name, data_dir = created

        Threads.@spawn run_ingest!(db, exp_id)

        log_action!(db, req; action = "experiment_created",
            entity_type = "experiment", entity_id = exp_id)

        HTTP.Response(202, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => exp_id, :status => "scanning",
                             :name => exp_name, :data_dir => data_dir)))
    end

    @get "/api/experiments" function(req::HTTP.Request)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM experiments ORDER BY id"))
        out = [_experiment_row_to_json(r, db) for r in rows]
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
    end

    @get "/api/experiments/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(_experiment_row_to_json(rows[1], db)))
    end

    @patch "/api/experiments/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))

        body = json(req)

        for k in _READONLY_FIELDS
            (haskey(body, Symbol(k)) || haskey(body, k)) &&
                return HTTP.Response(400,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "$k is read-only; change it via create/scan")))
        end
        unknown = sort!([String(k) for k in keys(body) if !(String(k) in _PLAIN_PATCH_FIELDS)])
        isempty(unknown) || return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "not an experiment field: $(join(unknown, ", ")); " *
                                       "only name and description can be changed")))

        set_clauses = String[]
        params      = Any[]
        changed     = Dict{Symbol, Any}()

        for field in _PLAIN_PATCH_FIELDS
            (haskey(body, Symbol(field)) || haskey(body, field)) || continue
            val = get(body, Symbol(field), get(body, field, nothing))
            val === nothing && !(field in _CLEARABLE_PATCH_FIELDS) && continue
            if field == "name"
                (val isa AbstractString && !isempty(strip(val))) || return HTTP.Response(400,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "name must be a non-empty string")))
                val = String(strip(val))
            end
            push!(set_clauses, "$field = ?"); push!(params, something(val, missing))
            changed[Symbol(field)] = val
        end

        isempty(set_clauses) && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "no patchable fields; supply name or description")))

        push!(params, id)
        return with_idempotency(db, req) do
            DBInterface.execute(db,
                "UPDATE experiments SET $(join(set_clauses, ", ")) WHERE id = ?",
                params)
            result = apply_event!(InTransaction(), db, req;
                kind = "experiment_updated",
                entity_type = "experiment", entity_id = id,
                payload = changed)
            _enqueue_broadcast_from_result!(result, "experiment_updated", "experiment", id)
            updated_rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT * FROM experiments WHERE id = ?", [id]))
            HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(_experiment_row_to_json(updated_rows[1], db)))
        end
    end

    @get "/api/experiments/{id}/problems" function(req::HTTP.Request, id::Int)
        db = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT ingest_error FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return _json(404, Dict(:error => "experiment not found"))
        _json(200, Dict(:ingest_error => coalesce(first(rows).ingest_error, nothing),
                        :exposures => ingest_problems(db, id),
                        :setups => [Dict(:root => r, :folder => basename(rstrip(r, '/')))
                                    for r in unrecognized_setup_roots(db, id)]))
    end

    @post "/api/experiments/{id}/setups" function(req::HTTP.Request, id::Int)
        db = current_db()
        body = json(req)
        _is_obj(body) || return _json_error(400, Dict(:error => "body must be an object"))
        added = _add_beamtime_setups(db, id, body)
        added isa HTTP.Response && return added
        _set_ingest_status!(db, id, "scanning")
        Threads.@spawn run_ingest!(db, id; phase = "rescan")
        log_action!(db, req; action = "setups_added", entity_type = "experiment", entity_id = id)
        HTTP.Response(202, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:status => "scanning", :experiment_id => id, :setup_ids => added)))
    end

    @post "/api/experiments/{id}/scan" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, ingest_status FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))

        req_body = isempty(req.body) ? Dict{String,Any}() : try
            JSON3.read(req.body)
        catch
            Dict{String,Any}()
        end
        force_scan = get(req_body, :force, get(req_body, "force", false)) === true
        Threads.@spawn with_background_db(bg -> try_backfill_monitors!(bg, id), db)
        if !(force_scan || ingest_running(id) || ingest_has_work(db, id))
            return HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:status => coalesce(first(rows).ingest_status, "complete"),
                                 :experiment_id => id, :changed => false)))
        end
        phase = _ingest_phase(db, id)
        _set_ingest_status!(db, id, "scanning")
        Threads.@spawn run_ingest!(db, id; phase, force = force_scan)

        log_action!(db, req; action = "scan",
            entity_type = "experiment", entity_id = id)

        HTTP.Response(202, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:status => "scanning", :experiment_id => id)))
    end

    @delete "/api/experiments/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, ingest_status FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))
        isequal(first(rows).ingest_status, "scanning") && return HTTP.Response(409,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "the experiment is still ingesting")))


        # Live schema (db.jl:32-152, confirmed 2026-06-18) keys a DEEP tree of
        # structural tables off exposures/samples/indices, and almost none of those
        # FKs declare `ON DELETE CASCADE`: samples.experiment_id (db.jl:34),
        # exposures.sample_id (db.jl:42), and the exposure-keyed tables
        # exposure_sources, indices, auto_peaks, peak_curations,
        # index_groups, assignments, assignment_members (db.jl:60-152), plus the
        # indices-keyed index_peaks / index_group_members. Corrected Plan A does NOT
        # add cascades (SQLite cannot ALTER a cascade onto an existing column).
        # FK enforcement is ON at the connection level (open_db, db.jl:1907), so a
        # bare `DELETE FROM experiments` would FK-fail and 500.
        #
        # Enumerating ~16 child deletes in FK order is brittle and easy to leave
        # incomplete as the schema grows. Instead follow the codebase's own teardown
        # idiom for cross-FK structural surgery (db.jl:1620-1647 and the other
        # migrations): toggle `PRAGMA foreign_keys = OFF` at the CONNECTION level
        # OUTSIDE the transaction (it is a documented no-op mid-transaction), do the
        # parent + cascade deletes, then restore `ON` in a `finally`. Delete in
        # FK-child→parent order anyway so the row set is internally consistent.
        #
        # Concurrency note: this disables FK checks connection-wide for the duration
        # of the delete on the shared singleton connection (parallel=true). Acceptable
        # because (a) experiment delete is a rare admin action, (b) the whole delete
        # is serialized under `_DB_WRITE_LOCK`, and (c) it mirrors the established
        # migration precedent. Do NOT issue `PRAGMA foreign_keys` INSIDE the
        # transaction — SQLite silently ignores it there.
        lock(_DB_WRITE_LOCK) do
            DBInterface.execute(db, "PRAGMA foreign_keys = OFF")
            try
                SQLite.transaction(db) do
                    # exposure-keyed structural tables (children of exposures/indices)
                    DBInterface.execute(db, """
                        DELETE FROM index_peaks WHERE index_id IN
                          (SELECT i.id FROM indices i
                             JOIN exposures e ON e.id = i.exposure_id
                            WHERE e.experiment_id = ?)""", [id])
                    DBInterface.execute(db, """
                        DELETE FROM index_group_members WHERE index_id IN
                          (SELECT i.id FROM indices i
                             JOIN exposures e ON e.id = i.exposure_id
                            WHERE e.experiment_id = ?)""", [id])
                    for tbl in ("assignment_members", "assignments", "index_groups",
                                "indices", "auto_peaks", "peak_curations",
                                "exposure_sources")
                        # exposure_sources keys two exposure columns; clean both.
                        if tbl == "exposure_sources"
                            DBInterface.execute(db, """
                                DELETE FROM exposure_sources WHERE averaged_exposure_id IN
                                  (SELECT id FROM exposures WHERE experiment_id = ?)
                                   OR source_exposure_id IN
                                  (SELECT id FROM exposures WHERE experiment_id = ?)""",
                                [id, id])
                        else
                            DBInterface.execute(db, """
                                DELETE FROM $tbl WHERE exposure_id IN
                                  (SELECT id FROM exposures WHERE experiment_id = ?)""",
                                [id])
                        end
                    end
                    # sample-keyed tables, then the core rows.
                    DBInterface.execute(db, """
                        DELETE FROM sample_bridges WHERE sample_id IN
                          (SELECT id FROM samples WHERE experiment_id = ?)""", [id])
                    DBInterface.execute(db, """
                        DELETE FROM sample_compares WHERE sample_id IN
                          (SELECT id FROM samples WHERE experiment_id = ?)""", [id])
                    DBInterface.execute(db, """
                        DELETE FROM sample_subtractions WHERE sample_id IN
                          (SELECT id FROM samples WHERE experiment_id = ?)""", [id])
                    DBInterface.execute(db, """
                        DELETE FROM sample_labels WHERE sample_id IN
                          (SELECT id FROM samples WHERE experiment_id = ?)""", [id])
                    DBInterface.execute(db, """
                        DELETE FROM sample_calls WHERE sample_id IN
                          (SELECT id FROM samples WHERE experiment_id = ?)""", [id])
                    DBInterface.execute(db,
                        "DELETE FROM exposures WHERE experiment_id = ?", [id])
                    DBInterface.execute(db,
                        "DELETE FROM samples WHERE experiment_id = ?", [id])
                    DBInterface.execute(db,
                        "DELETE FROM loads WHERE experiment_id = ?", [id])
                    DBInterface.execute(db,
                        "DELETE FROM setups WHERE experiment_id = ?", [id])
                    DBInterface.execute(db,
                        "DELETE FROM experiments WHERE id = ?", [id])
                end
            finally
                DBInterface.execute(db, "PRAGMA foreign_keys = ON")
            end
        end

        log_action!(db, req; action = "experiment_deleted",
            entity_type = "experiment", entity_id = id)

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => id, :deleted => true)))
    end
end
