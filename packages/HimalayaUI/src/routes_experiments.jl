using HTTP, JSON3, DBInterface, Tables, Oxygen

# Mutable geometry fields. Each writable field has a companion *_source column
# (set to "user" on override; never refreshed by rescan). name/description and
# the path fields (data_dir, analysis_dir, manifest_path) remain read-only here —
# name/description editing lands in Phase E1 (needs a description-column migration);
# path fields are set at create time and must stay in sync with the filesystem.
const _GEOMETRY_PATCH_FIELDS = [
    "flight_path_m", "beam_center_x", "beam_center_y",
    "pixel_size_um", "energy_kev", "q_units",
]
const _READONLY_FIELDS = ["data_dir", "analysis_dir", "manifest_path", "path",
                           "id", "created_at", "name", "description"]

"""
    _beamline_from_config(cfg_text) -> NamedTuple

Extract the render-only `[beamline]` fields from an experiment's stored config
TOML: `q_units` (String, default "A-1") and `beam_center_x`, `beam_center_y`,
`pixel_size_um` (each `Float64` or `nothing`). A bare integer in the TOML is
coerced to `Float64` to match the REAL-column behaviour of energy_kev.

Defensive by design: a non-string, empty, malformed, or non-numeric value all
fall back to the all-default tuple — one experiment's bad config must never 500
a list endpoint. The numeric coercion lives INSIDE the try so a value like
`beam_center_x = "oops"` throws into the fallback rather than out of the route.
"""
function _beamline_from_config(cfg_text)
    default = (q_units = "A-1", beam_center_x = nothing,
               beam_center_y = nothing, pixel_size_um = nothing)
    if cfg_text isa AbstractString && !isempty(cfg_text)
        try
            bl = get(TOML.parse(cfg_text), "beamline", Dict())
            num(k) = (v = get(bl, k, nothing); v === nothing ? nothing : Float64(v))
            # Guard the q_units String contract: a non-string TOML value (e.g.
            # `q_units = 5`) would otherwise satisfy the NamedTuple but throw at
            # the shim's `::String` annotation — OUTSIDE this try — and 500 the
            # samples list route. Coerce to the default here so the docstring's
            # "never 500 a list endpoint" holds for q_units too.
            qu = get(bl, "q_units", "A-1")
            return (q_units       = qu isa AbstractString ? qu : "A-1",
                    beam_center_x = num("beam_center_x"),
                    beam_center_y = num("beam_center_y"),
                    pixel_size_um = num("pixel_size_um"))
        catch
            return default
        end
    end
    return default
end

# Back-compat shim: routes_samples.jl still calls this for its per-sample q_units.
_q_units_from_config(cfg_text)::String = _beamline_from_config(cfg_text).q_units

"""
    _experiment_stats(db, exp_id) -> NamedTuple

Cheap roll-up of counts for the shared experiment header stat ledger.
"""
function _experiment_stats(db::SQLite.DB, exp_id::Integer)
    loads = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM loads WHERE experiment_id = ?", [exp_id]))[1].c
    samples = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM samples WHERE experiment_id = ?", [exp_id]))[1].c
    exposures = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM exposures WHERE experiment_id = ?", [exp_id]))[1].c
    (loads = Int(loads), samples = Int(samples), exposures = Int(exposures))
end

"""
    _experiment_row_to_json(row, db) -> Dict

Serialize an experiments row to the wire format. Reads typed geometry columns
from Phase A directly (no TOML overlay). Falls back to `_beamline_from_config`
for legacy rows that still have their geometry only in the TOML `config` blob
(experiments ingested before Phase A).
"""
function _experiment_row_to_json(row::NamedTuple, db::Union{SQLite.DB, Nothing} = nothing)
    d = row_to_json(row)
    # Prefer typed columns (Phase A); fall back to TOML blob for legacy rows.
    has_typed = !isnothing(get(d, :beam_center_x, nothing)) ||
                !isnothing(get(d, :energy_kev, nothing))
    if !has_typed
        bl = _beamline_from_config(get(d, :config, nothing))
        d[:q_units]              = bl.q_units
        d[:beam_center_x]        = bl.beam_center_x
        d[:beam_center_y]        = bl.beam_center_y
        d[:pixel_size_um]        = bl.pixel_size_um
        # energy_kev / flight_path_m are real columns even pre-Phase-A
        # (live create_experiment! writes them); surface their VALUE keys too so
        # the wire shape is identical to the typed path (a legacy row simply
        # reports whatever those columns hold — possibly nothing — never absent).
        d[:energy_kev]           = get(d, :energy_kev, nothing)
        d[:flight_path_m]        = get(d, :flight_path_m, nothing)
        d[:beam_center_x_source] = "default"
        d[:beam_center_y_source] = "default"
        d[:pixel_size_um_source] = "default"
        d[:energy_kev_source]    = "default"
        d[:flight_path_m_source] = "default"
        d[:q_units_source]       = "default"
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

        # Required: path to the data directory.
        path_val = get(body, :path, get(body, "path", nothing))
        path_val === nothing && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "path is required")))
        data_dir = String(path_val)

        isdir(data_dir) || return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "path does not exist or is not a directory",
                             :path  => data_dir)))

        # Derive defaults.
        name_val  = get(body, :name, get(body, "name", nothing))
        exp_name  = name_val !== nothing ? String(name_val) : basename(rstrip(data_dir, '/'))
        # analysis_dir convention: look for an `analysis` subdirectory; fall back to data_dir.
        analysis_dir = let ad = joinpath(data_dir, "analysis")
            isdir(ad) ? ad : data_dir
        end

        exp_id = lock(_DB_WRITE_LOCK) do
            SQLite.transaction(db) do
                create_experiment!(db;
                    name         = exp_name,
                    path         = data_dir,
                    data_dir     = data_dir,
                    analysis_dir = analysis_dir,
                    ingest_status = "scanning")
            end
        end

        broadcast_progress!(exp_id; kind = "ingest_started", processed = 0, total = 0)

        # Kick off first scan asynchronously.
        Threads.@spawn begin
            try
                # scan_and_group! (Phase B, ingest.jl) resolves data_dir from the row
                # and is idempotent (dedup INSERT keys), so first-scan == rescan.
                scan_and_group!(db, exp_id)
                lock(_DB_WRITE_LOCK) do
                    SQLite.transaction(db) do
                        DBInterface.execute(db,
                            "UPDATE experiments SET ingest_status = 'complete', last_scanned_at = ? WHERE id = ?",
                            [format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ"), exp_id])
                    end
                end
                broadcast_progress!(exp_id; kind = "ingest_complete", processed = 0, total = 0)
                # Arm the rescan scheduler after a successful first scan.
                start_rescan_scheduler!(db, exp_id)
            catch err
                @warn "first scan failed" experiment_id = exp_id exception = err
                lock(_DB_WRITE_LOCK) do
                    SQLite.transaction(db) do
                        DBInterface.execute(db,
                            "UPDATE experiments SET ingest_status = 'failed' WHERE id = ?", [exp_id])
                    end
                end
                broadcast_progress!(exp_id; kind = "ingest_failed",
                    processed = 0, total = 0, error = sprint(showerror, err))
            end
        end

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
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write([_experiment_row_to_json(r) for r in rows]))
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

        # No recognized patchable (geometry) field in the body → 400, matching the
        # codebase's PATCH validation convention (cf. PATCH /samples, exercised by
        # test_route_validation_routing.jl: a body with no patchable field is a
        # bad request, not a 200 no-op).
        isempty(set_clauses) && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "no patchable fields; supply a geometry field (flight_path_m, beam_center_x/y, pixel_size_um, energy_kev, q_units)")))

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

        loads_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM loads WHERE experiment_id = ? ORDER BY load_index", [id]))

        result = map(loads_rows) do lr
            samples_rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT * FROM samples WHERE load_id = ? ORDER BY slot_index", [Int(lr.id)]))

            samples = map(samples_rows) do sr
                exposures_rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, filename, prp_path, timestamp, exposure_time,
                            horizontal_position, scan_id, frame_no, status, selected,
                            image_path, content_fingerprint
                       FROM exposures WHERE sample_id = ? ORDER BY frame_no, id",
                    [Int(sr.id)]))
                d = row_to_json(sr)
                # exposures.selected is a SQLite 0/1 int; coerce to a JSON bool so this
                # endpoint matches every other exposure serializer (routes_analysis.jl /
                # comparisons.jl ~493 use `row_to_json(e; bool_keys=(:selected,))`).
                d[:exposures] = [row_to_json(e; bool_keys = (:selected,)) for e in exposures_rows]
                d
            end

            d = row_to_json(lr)
            d[:samples] = samples
            d
        end

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(result))
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
            "SELECT id FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))

        # Mark as scanning immediately so the frontend header can show progress.
        lock(_DB_WRITE_LOCK) do
            SQLite.transaction(db) do
                DBInterface.execute(db,
                    "UPDATE experiments SET ingest_status = 'scanning' WHERE id = ?", [id])
            end
        end

        broadcast_progress!(id; kind = "ingest_started", processed = 0, total = 0)

        # Run the cheap change-check + additive scan on a @spawn'd task so this request
        # returns immediately; progress streams over SSE. Both Phase B functions
        # (cheap_change_check, scan_and_group!) resolve the experiment's data_dir from
        # the row themselves, and return gracefully on an empty directory.
        Threads.@spawn begin
            try
                changed = cheap_change_check(db, id)
                if changed
                    scan_and_group!(db, id)
                    start_rescan_scheduler!(db, id)   # re-arm the fast-tier scheduler
                end

                lock(_DB_WRITE_LOCK) do
                    SQLite.transaction(db) do
                        DBInterface.execute(db,
                            "UPDATE experiments SET ingest_status = 'complete', last_scanned_at = ? WHERE id = ?",
                            [format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ"), id])
                    end
                end
                broadcast_progress!(id; kind = "ingest_complete",
                    processed = 0, total = 0, changed = changed)
            catch err
                @warn "scan failed" experiment_id = id exception = err
                lock(_DB_WRITE_LOCK) do
                    SQLite.transaction(db) do
                        DBInterface.execute(db,
                            "UPDATE experiments SET ingest_status = 'failed' WHERE id = ?", [id])
                    end
                end
                broadcast_progress!(id; kind = "ingest_failed",
                    processed = 0, total = 0,
                    error = sprint(showerror, err))
            end
        end

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

        # Stop the rescan timer BEFORE the DB delete so the timer callback cannot
        # fire against a non-existent row. stop_rescan_scheduler! is defined in server.jl
        # and is a no-op when no timer is running for this id.
        stop_rescan_scheduler!(id)

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
                        DELETE FROM sample_messages WHERE sample_id IN
                          (SELECT id FROM samples WHERE experiment_id = ?)""", [id])
                    # Cross-feature tables OUTSIDE the core tree (db.jl:700-1060) that
                    # reference this experiment's samples/exposures. FK enforcement is
                    # OFF here, so their declared ON DELETE actions do NOT fire —
                    # replicate each by hand so no orphan rows remain:
                    #   series_samples.sample_id   -> samples   ON DELETE CASCADE  (delete row)
                    #   series_members.exposure_id -> exposures ON DELETE SET NULL (null the ref)
                    #   comparison_members.exposure_id -> exposures ON DELETE SET NULL (null the ref)
                    # Must run BEFORE the samples/exposures deletes below so the
                    # IN-subqueries still resolve the ids.
                    DBInterface.execute(db, """
                        DELETE FROM series_samples WHERE sample_id IN
                          (SELECT id FROM samples WHERE experiment_id = ?)""", [id])
                    DBInterface.execute(db, """
                        UPDATE series_members SET exposure_id = NULL WHERE exposure_id IN
                          (SELECT id FROM exposures WHERE experiment_id = ?)""", [id])
                    DBInterface.execute(db, """
                        UPDATE comparison_members SET exposure_id = NULL WHERE exposure_id IN
                          (SELECT id FROM exposures WHERE experiment_id = ?)""", [id])
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
