using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_trace_routes!()
    @get "/api/exposures/{id}/trace" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT e.filename, e.kind, x.analysis_dir
             FROM exposures e JOIN samples s ON s.id = e.sample_id
             JOIN experiments x ON x.id = s.experiment_id
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

        path = joinpath(String(row.analysis_dir), String(row.filename) * ".dat")
        isfile(path) || return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => ".dat file not found: $path")))

        q, I, σ = load_dat(path)
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:q => q, :I => I, :sigma => σ)))
    end
end
