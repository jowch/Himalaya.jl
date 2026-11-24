using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_trace_routes!()
    @get "/api/exposures/{id}/trace" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = sqlrows(DBInterface.execute(db,
            "SELECT e.filename, e.kind, e.transmitted_counts, e.exposure_time, e.counting_time,
                    e.trace_median, e.analysis_error, e.incident_monitor,
                    e.setup_id, st.analysis_dir, st.integration_pattern, st.trace_normalisation, st.mode
             FROM exposures e JOIN setups st ON st.id = e.setup_id
             WHERE e.id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "image not found")))
        row = rows[1]
        is_combined_kind(row.kind) && return _combined_trace_response(db, id)
        String(row.kind) == "file" || return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "trace not available for derived images")))
        row.filename === missing && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "image has no filename")))

        pattern = integration_pattern_for(row)
        filename = String(row.filename)
        path = resolve_trace_path(String(row.analysis_dir), pattern, filename)
        path === nothing && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => ".dat file not found for $(row.filename)")))

        q, I, σ = load_dat(path)
        monitors = Dict(:transmitted_counts => coalesce(row.transmitted_counts, nothing),
                        :exposure_time      => coalesce(row.exposure_time, nothing),
                        :counting_time      => coalesce(row.counting_time, nothing),
                        :normalised_by      => coalesce(row.trace_normalisation, NORMALISATION_UNKNOWN),
                        :beam_off           => is_beam_off(row.trace_median, row.analysis_error, row.incident_monitor),
                        :i2_collapsed       => is_i2_collapsed(row.transmitted_counts, row.incident_monitor, row.counting_time,
                                                               setup_transmission_median(db, row.setup_id)),
                        :mode               => coalesce(row.mode, nothing))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:q => q, :I => I, :sigma => σ, :monitors => monitors)))
    end
end

function _combined_trace_response(db::SQLite.DB, id::Int)
    q, I, σ = try
        load_combined_trace(db, id)
    catch err
        err isa ArgumentError || rethrow()
        return HTTP.Response(404, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => err.msg)))
    end
    HTTP.Response(200, ["Content-Type" => "application/json"],
        JSON3.write(Dict(:q => q, :I => I, :sigma => σ, :monitors => combined_monitors(db, id),
                         :sources => combined_source_ids(db, id))))
end
