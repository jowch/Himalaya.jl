using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_trace_routes!()
    @get "/api/exposures/{id}/trace" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT e.filename, e.kind, e.transmitted_counts, e.exposure_time, e.counting_time,
                    st.analysis_dir, st.integration_pattern
             FROM exposures e JOIN setups st ON st.id = e.setup_id
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

        pattern = integration_pattern_for(row)
        filename = String(row.filename)
        path = resolve_trace_path(String(row.analysis_dir), pattern, filename)
        path === nothing && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => ".dat file not found for $(row.filename)")))

        q, I, σ = load_dat(path)
        monitors = Dict(:transmitted_counts => coalesce(row.transmitted_counts, nothing),
                        :exposure_time      => coalesce(row.exposure_time, nothing),
                        :counting_time      => coalesce(row.counting_time, nothing))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:q => q, :I => I, :sigma => σ, :monitors => monitors)))
    end
end
