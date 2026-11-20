using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_trace_routes!()
    @get "/api/exposures/{id}/trace" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT e.filename, e.kind, x.id AS experiment_id, x.analysis_dir
             FROM exposures e JOIN experiments x ON x.id = e.experiment_id
             WHERE e.id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure not found")))
        row = rows[1]
        String(row.kind) == "file" || return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "trace not available for derived exposures")))
        row.filename === missing && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "exposure has no filename")))

        # Resolve via the shared resolver (exact per-frame name, else the per-
        # acquisition `_tot.dat` — migrated filenames are the full per-frame stem
        # while the trace total is named by the acquisition stem). See
        # resolve_trace_path.
        pattern = integration_pattern_for(db, Int(row.experiment_id))
        filename = String(row.filename)
        path = resolve_trace_path(String(row.analysis_dir), pattern, filename)
        path === nothing && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => ".dat file not found for $(row.filename)")))

        q, I, σ = load_dat(path)
        per_frame = basename(path) == basename(replace(pattern, "{name}" => filename))
        frames = trace_frames(db, Int(row.experiment_id), id, filename, per_frame)
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:q => q, :I => I, :sigma => σ, :frames => frames)))
    end
end

function trace_frames(db::SQLite.DB, experiment_id::Int, exposure_id::Int,
                      filename::AbstractString, per_frame::Bool)
    cols = "id, transmitted_counts, exposure_time, counting_time"
    rows = if per_frame
        Tables.rowtable(DBInterface.execute(db,
            "SELECT $cols FROM exposures WHERE id = ?", [exposure_id]))
    else
        acquisition = acquisition_stem(filename)
        filter(r -> acquisition_stem(String(r.filename)) == acquisition,
            Tables.rowtable(DBInterface.execute(db,
                "SELECT $cols, filename FROM exposures
                  WHERE experiment_id = ? AND kind = 'file' AND substr(filename, 1, ?) = ?
                  ORDER BY id",
                [experiment_id, length(acquisition), acquisition])))
    end
    [Dict(:exposure_id => Int(r.id),
          :transmitted_counts => coalesce(r.transmitted_counts, nothing),
          :exposure_time => coalesce(r.exposure_time, nothing),
          :counting_time => coalesce(r.counting_time, nothing)) for r in rows]
end
