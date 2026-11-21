using HTTP, JSON3, DBInterface, Tables, Oxygen

# PATCH /experiments writes three tiers of fields into one UPDATE:
#  - geometry (_GEOMETRY_PATCH_FIELDS): each has a companion *_source column set
#    to "user" on override and never refreshed by rescan.
#  - plain (_PLAIN_PATCH_FIELDS = name, description): plain write, no *_source, no rescan.
#  - patterns (_PATTERN_FIELDS): plain write; the client forces a rescan afterwards
#    (the discovered file set changes when a pattern changes).
# The path fields (data_dir, analysis_dir, manifest_path, path) + id/created_at stay
# read-only (_READONLY_FIELDS): set at create time, must stay in sync with the filesystem.
const _GEOMETRY_PATCH_FIELDS = [
    "flight_path_m", "beam_center_x", "beam_center_y",
    "pixel_size_um", "energy_kev", "q_units",
]
const _READONLY_FIELDS = ["data_dir", "analysis_dir", "manifest_path", "path",
                           "id", "created_at"]
const _PLAIN_PATCH_FIELDS = ["name", "description"]
const _CLEARABLE_PATCH_FIELDS = ("description",)
const _PATTERN_FIELDS     = ["image_pattern", "metadata_pattern", "integration_pattern"]

# Serialize a GroupingFlag (MergeFlag/SplitFlag) to the §8.8 wire shape. The
# internal structs carry no `kind` field, but the frontend `GroupingFlag` union
# narrows on `kind` — inject it here at the wire boundary (GET /experiments/:id/loads).
_flag_json(::Nothing) = nothing
_flag_json(f::MergeFlag) = Dict(:kind => "merge",
    :merge_with_sample_id => f.merge_with_sample_id, :merge_with_label => f.merge_with_label)
_flag_json(f::SplitFlag) = Dict(:kind => "split",
    :split_at_index => f.split_at_index, :jump_from => f.jump_from, :jump_to => f.jump_to)

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

"""
    _experiment_review_count(db, exp_id) -> Int

Count of samples in this experiment that carry a pending grouping flag
(SplitFlag or MergeFlag) — i.e. samples the user should review. Reuses the
single-sourced flag derivation in `get_loads_rollup` / `derive_sample_flags`.

# ponytail: O(rollup) per experiment, fine for the gallery (loads once, few
# experiments); precompute/lighter-query if it ever slows at many-experiments scale.
"""
function _experiment_review_count(db::SQLite.DB, exp_id::Integer)
    rollup = get_loads_rollup(db, exp_id)   # nested loads → samples (each .flag); dismiss-suppressed
    count(s -> s.flag !== nothing, (s for ld in rollup for s in ld.samples))
end

"""
    _experiment_row_to_json(row, db) -> Dict

Serialize an experiments row to the wire format.
"""
function _experiment_row_to_json(row::NamedTuple, db::Union{SQLite.DB, Nothing} = nothing)
    d = row_to_json(row)
    delete!(d, :config)
    # q_units is a unit convention, never derived from the data files (geometry.jl
    # `derive_geometry` produces no q_units). Default it to the SAXS convention
    # "A-1" whenever the column is unset — including the typed path above, where a
    # migrated/fresh row carries real geometry but a NULL q_units — so every
    # experiment reports a concrete unit (the geometry ledger never shows it blank)
    # and stays consistent with routes_samples' per-sample A-1 default. A user
    # PATCH override leaves the column non-null, so this is skipped.
    let qv = get(d, :q_units, nothing)
        if isnothing(qv) || ismissing(qv)
            d[:q_units]        = "A-1"
            d[:q_units_source] = "default"
        end
    end
    # Add stats roll-up when db is supplied (single-row endpoint).
    if db !== nothing
        exp_id = Int(row.id)
        d[:stats] = _experiment_stats(db, exp_id)
    end
    d
end

function register_experiments_routes!()
    @post "/api/experiments" function(req::HTTP.Request)
        db   = current_db()
        body = json(req)

        # Required: the data directory. The funnel resolver (/api/fs/resolve) sends
        # an explicit `data_dir` (+ `analysis_dir`, `name`) that the user confirmed
        # in Configuration; older callers send `path`. Prefer the explicit field.
        path_val = get(body, :data_dir, get(body, "data_dir",
                       get(body, :path, get(body, "path", nothing))))
        path_val === nothing && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "data_dir (or path) is required")))
        data_dir = String(path_val)

        isdir(data_dir) || return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "path does not exist or is not a directory",
                             :path  => data_dir)))

        # One experiment per directory: an experiment IS its data_dir (ingest
        # everything in it), so a second experiment on the same dir would be a
        # duplicate. Reject (the UI also validates this up front).
        existing_dir = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM experiments WHERE data_dir = ?", [data_dir]))
        isempty(existing_dir) || return HTTP.Response(409,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "an experiment already uses this directory",
                             :path  => data_dir,
                             :experiment_id => Int(first(existing_dir).id))))

        # Derive defaults.
        name_val  = get(body, :name, get(body, "name", nothing))
        exp_name  = name_val !== nothing ? String(name_val) : basename(rstrip(data_dir, '/'))
        # analysis_dir: prefer the explicit, user-confirmed value from the funnel
        # (the resolver's default that the user corrected); else the old convention
        # (look for an `analysis` subdir, fall back to data_dir) for older callers.
        ad_val = get(body, :analysis_dir, get(body, "analysis_dir", nothing))
        analysis_dir = ad_val !== nothing ? String(ad_val) :
            let ad = joinpath(data_dir, "analysis"); isdir(ad) ? ad : data_dir end

        # Patterns are edited on Configuration before Approve; persisted to the row
        # so scan_and_group! (which reads the row + coalesces) uses them. Nested
        # under `patterns` to match the frontend CreateExperimentBody.
        pats = get(body, :patterns, get(body, "patterns", Dict()))
        ppat(k) = (v = get(pats, k, get(pats, string(k), nothing)); v === nothing ? nothing : String(v))

        # Geometry (optional): the funnel preview already derived geometry with the
        # confirmed setup file, so Approve sends the WHOLE geometry — each field as
        # a value + its honest source (`setup`/`prp`/`computed`, or `user` for a
        # manual edit). It is persisted verbatim and the scan's fill-only step
        # (`_fill_unset_geometry!`, ingest.jl) leaves any established field alone, so
        # geometry is derived once (at preview) and never re-derived. A field with a
        # value but no source defaults to `user` (back-compat with a values-only
        # body). Non-numeric values / non-string sources are ignored (treated as not
        # provided) so a malformed body can't 500 the create.
        geo_body = get(body, :geometry, get(body, "geometry", Dict()))
        gnum(k)  = (v = get(geo_body, k, get(geo_body, string(k), nothing)); v isa Number ? Float64(v) : nothing)
        gstr(k)  = (v = get(geo_body, k, get(geo_body, string(k), nothing)); v isa AbstractString ? String(v) : nothing)
        # value present → its source (or 'user' when omitted); absent → 'default'.
        psrc(val, src) = val === nothing ? "default" : (src === nothing ? "user" : src)
        g_bcx = gnum(:beam_center_x); g_bcy = gnum(:beam_center_y)
        g_fp  = gnum(:flight_path_m); g_px  = gnum(:pixel_size_um); g_en = gnum(:energy_kev)

        exp_id = lock(_DB_WRITE_LOCK) do
            SQLite.transaction(db) do
                create_experiment!(db;
                    name              = exp_name,
                    path              = data_dir,
                    data_dir          = data_dir,
                    analysis_dir      = analysis_dir,
                    image_pattern     = ppat(:image),
                    metadata_pattern  = ppat(:metadata),
                    integration_pattern = ppat(:integration),
                    beam_center_x = g_bcx, beam_center_x_source = psrc(g_bcx, gstr(:beam_center_x_source)),
                    beam_center_y = g_bcy, beam_center_y_source = psrc(g_bcy, gstr(:beam_center_y_source)),
                    flight_path_m = g_fp,  flight_path_m_source = psrc(g_fp,  gstr(:flight_path_m_source)),
                    pixel_size_um = g_px,  pixel_size_um_source = psrc(g_px,  gstr(:pixel_size_um_source)),
                    energy_kev    = g_en,  energy_kev_source    = psrc(g_en,  gstr(:energy_kev_source)),
                    ingest_status     = "scanning")
            end
        end

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
        out = map(rows) do r
            d = _experiment_row_to_json(r, db)
            d[:review_count] = _experiment_review_count(db, Int(r.id))
            d
        end
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

        # Reject any attempt to write read-only path/id fields.
        for k in _READONLY_FIELDS
            (haskey(body, Symbol(k)) || haskey(body, k)) &&
                return HTTP.Response(400,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "$k is read-only; change it via create/scan")))
        end

        # Build SET clauses for geometry fields present in the body.
        set_clauses = String[]
        params      = Any[]
        for field in _GEOMETRY_PATCH_FIELDS
            val = get(body, Symbol(field), get(body, field, nothing))
            val === nothing && continue
            push!(set_clauses, "$field = ?")
            push!(params, val)
            push!(set_clauses, "$(field)_source = 'user'")
        end

        for field in _PLAIN_PATCH_FIELDS
            (haskey(body, Symbol(field)) || haskey(body, field)) || continue
            val = get(body, Symbol(field), get(body, field, nothing))
            val === nothing && !(field in _CLEARABLE_PATCH_FIELDS) && continue
            push!(set_clauses, "$field = ?"); push!(params, something(val, missing))
        end

        for field in _PATTERN_FIELDS
            val = get(body, Symbol(field), get(body, field, nothing))
            val === nothing && continue
            push!(set_clauses, "$field = ?"); push!(params, val)
        end

        # No recognized patchable field in the body → 400, matching the
        # codebase's PATCH validation convention (cf. PATCH /samples, exercised by
        # test_route_validation_routing.jl: a body with no patchable field is a
        # bad request, not a 200 no-op).
        isempty(set_clauses) && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "no patchable fields; supply name, description, a pattern field (image_pattern, metadata_pattern, integration_pattern), or a geometry field (flight_path_m, beam_center_x/y, pixel_size_um, energy_kev, q_units)")))

        push!(params, id)
        lock(_DB_WRITE_LOCK) do
            SQLite.transaction(db) do
                DBInterface.execute(db,
                    "UPDATE experiments SET $(join(set_clauses, ", ")) WHERE id = ?",
                    params)
            end
        end

        updated_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM experiments WHERE id = ?", [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(_experiment_row_to_json(updated_rows[1], db)))
    end

    @get "/api/experiments/{id}/loads" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))

        rollup = get_loads_rollup(db, id)
        # Serialize: NamedTuples → JSON-friendly Dict tree.
        # Field names are the §8.8 contract verbatim (Load / LoadSample /
        # LoadExposure). `flag` is a GroupingFlag struct (MergeFlag/SplitFlag) or
        # nothing; `_flag_json` injects the `kind` discriminator the frontend
        # `GroupingFlag` union narrows on (the bare struct has no `kind` field),
        # emitting {"kind":"merge",…}/{"kind":"split",…} or JSON null.
        out = map(rollup) do ld
            Dict(
                :load_id     => ld.load_id,
                :load_index  => ld.load_index,
                :session_id  => ld.session_id,
                :start_time  => ld.start_time,
                :end_time    => ld.end_time,
                :frame_count => ld.frame_count,
                :note        => ld.note,
                :samples     => map(ld.samples) do sm
                    Dict(
                        :sample_id       => sm.sample_id,
                        :name            => sm.name,
                        :slot_index      => sm.slot_index,
                        :grouping_source => sm.grouping_source,
                        :name_source     => sm.name_source,
                        :merged_into_id  => sm.merged_into_id,
                        :flag            => _flag_json(sm.flag),
                        :exposures       => map(sm.exposures) do ex
                            Dict(
                                :id                  => ex.id,
                                :filename            => ex.filename,
                                :horizontal_position => ex.horizontal_position,
                                :timestamp           => ex.timestamp)
                        end)
                end)
        end
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
    end

    @post "/api/experiments/{id}/analyze" function(req::HTTP.Request, id::Int)
        db   = current_db()
        # Verify the experiment exists before iterating its samples.
        exists = !isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM experiments WHERE id = ?", [id])))
        exists || return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))

        samples   = get_samples(db, id)
        analyzed  = 0
        skipped   = String[]
        for sm in samples
            for ex in get_exposures(db, Int(sm.id))
                try
                    analyze_exposure!(db, Int(ex.id))
                    analyzed += 1
                catch e
                    push!(skipped, "$(sm.name)/$(ex.filename): $(sprint(showerror, e))")
                end
            end
        end

        log_action!(db, req; action = "analyze",
            entity_type = "experiment", entity_id = id,
            note = "analyzed=$analyzed skipped=$(length(skipped))")

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:analyzed => analyzed, :skipped => skipped)))
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
            "SELECT id FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))


        # Live schema (db.jl:32-152, confirmed 2026-06-18) keys a DEEP tree of
        # structural tables off exposures/samples/indices, and almost none of those
        # FKs declare `ON DELETE CASCADE`: samples.experiment_id (db.jl:34),
        # exposures.sample_id (db.jl:42), and the exposure-keyed tables
        # exposure_sources, exposure_tags, indices, auto_peaks, peak_curations,
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
                                "exposure_sources", "exposure_tags")
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
                        DELETE FROM sample_tags WHERE sample_id IN
                          (SELECT id FROM samples WHERE experiment_id = ?)""", [id])
                    DBInterface.execute(db, """
                        DELETE FROM sample_bridges WHERE sample_id IN
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
