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

        isempty(set_clauses) && return HTTP.Response(200,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => id, :updated => false)))

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
end
