using HTTP, JSON3, DBInterface, Tables, Oxygen

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

function _experiment_row_to_json(row::NamedTuple)
    d  = row_to_json(row)
    bl = _beamline_from_config(get(d, :config, nothing))
    d[:q_units]       = bl.q_units
    d[:beam_center_x] = bl.beam_center_x
    d[:beam_center_y] = bl.beam_center_y
    d[:pixel_size_um] = bl.pixel_size_um
    d
end

function register_experiments_routes!()
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
            JSON3.write(_experiment_row_to_json(rows[1])))
    end

    @patch "/api/experiments/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)

        # Experiment name and path fields are no longer mutable via PATCH.
        # Name is derived from experiment.toml via reingest; path fields
        # (data_dir, analysis_dir, manifest_path) must also go through reingest
        # to stay in sync with the config blob.
        # This route is a defensive surface for future fields only.
        return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment metadata is read-only; rename via experiment.toml + reingest")))
    end

    @post "/api/experiments/{id}/analyze" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT analysis_dir FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))
        analysis_dir = rows[1].analysis_dir

        samples   = get_samples(db, id)
        analyzed  = 0
        skipped   = String[]
        for sm in samples
            for ex in get_exposures(db, Int(sm.id))
                try
                    analyze_exposure!(db, Int(ex.id), String(analysis_dir))
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

    @post "/api/experiments/{id}/reingest" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT path FROM experiments WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "experiment not found")))
        exp_path = String(rows[1].path)
        try
            res = reingest!(db, id, exp_path)
            log_action!(db, req; action = "reingest",
                entity_type = "experiment", entity_id = id)
            return HTTP.Response(200,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:status          => String(res.status),
                                 :added_samples   => res.added_samples,
                                 :added_exposures => res.added_exposures,
                                 :manifest_path   => res.manifest_path)))
        catch e
            if e isa ManifestValidationError
                return HTTP.Response(400,
                    ["Content-Type" => "application/json"],
                    JSON3.write(Dict(:error => "manifest_invalid",
                                     :violations => [Dict(:kind => string(v.kind),
                                                          :sample_index => v.sample_index,
                                                          :sample_name => v.sample_name,
                                                          :detail => v.detail)
                                                     for v in e.violations])))
            end
            return HTTP.Response(500,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => sprint(showerror, e))))
        end
    end
end
