using HTTP, JSON3, DBInterface, Tables, Oxygen, CSV, SQLite

"""
    build_export(db, experiment_id) -> Vector{Dict}

Compose the per-sample summary used by both JSON and CSV exports.
For each exposure, `indices` reflects the active group if one is set,
otherwise all candidate indices ordered by score DESC.
"""
function build_export(db::SQLite.DB, experiment_id::Int)
    samples = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM samples WHERE experiment_id = ? ORDER BY id",
        [experiment_id]))

    out = Dict[]
    for sm in samples
        sid  = Int(sm.id)
        tags = Tables.rowtable(DBInterface.execute(db,
            "SELECT key, value FROM sample_tags WHERE sample_id = ?", [sid]))

        exposures = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM exposures WHERE sample_id = ? ORDER BY id", [sid]))

        ex_out = Dict[]
        for e in exposures
            eid = Int(e.id)
            active = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, kind FROM index_groups
                 WHERE exposure_id = ? AND active = 1
                 ORDER BY id DESC LIMIT 1", [eid]))
            ag_kind   = isempty(active) ? nothing : String(active[1].kind)
            active_id = isempty(active) ? nothing : Int(active[1].id)

            indices = if active_id === nothing
                Tables.rowtable(DBInterface.execute(db,
                    "SELECT phase, basis, score, r_squared, lattice_d
                     FROM indices WHERE exposure_id = ? ORDER BY score DESC",
                    [eid]))
            else
                Tables.rowtable(DBInterface.execute(db,
                    "SELECT i.phase, i.basis, i.score, i.r_squared, i.lattice_d
                     FROM indices i
                     JOIN index_group_members m ON m.index_id = i.id
                     WHERE m.group_id = ? ORDER BY i.score DESC", [active_id]))
            end

            push!(ex_out, Dict(
                :filename          => e.filename,
                :selected          => e.selected != 0,
                :active_group_kind => ag_kind,
                :indices           => rows_to_json(indices),
            ))
        end

        push!(out, Dict(
            :label     => sm.label,
            :name      => sm.name,
            :notes     => sm.notes isa Missing ? nothing : sm.notes,
            :tags      => rows_to_json(tags),
            :exposures => ex_out,
        ))
    end
    out
end

function _export_csv(summary)
    header = ["sample_label","sample_name","exposure_filename",
              "phases","lattice_ds","r_squareds","scores","tags"]
    rows = Vector{Vector{Any}}()
    for sm in summary
        tag_str = join(["$(t[:key])=$(t[:value])" for t in sm[:tags]], ";")
        if isempty(sm[:exposures])
            push!(rows, [sm[:label], sm[:name], "", "", "", "", "", tag_str])
            continue
        end
        for ex in sm[:exposures]
            phases = join([String(ix[:phase])     for ix in ex[:indices]], ";")
            lds    = join([string(ix[:lattice_d]) for ix in ex[:indices]], ";")
            r2s    = join([string(ix[:r_squared]) for ix in ex[:indices]], ";")
            scs    = join([string(ix[:score])     for ix in ex[:indices]], ";")
            push!(rows, [sm[:label], sm[:name], ex[:filename],
                         phases, lds, r2s, scs, tag_str])
        end
    end
    io = IOBuffer()
    CSV.write(io, (; (Symbol(h) => [r[i] for r in rows] for (i, h) in enumerate(header))...))
    String(take!(io))
end

function register_export_routes!()
    @get "/api/experiments/{id}/export" function(req::HTTP.Request, id::Int)
        fmt = try
            q = HTTP.queryparams(HTTP.URI(req.target))
            get(q, "format", "json")
        catch
            "json"
        end

        db      = current_db()
        summary = build_export(db, id)

        if fmt == "json"
            return HTTP.Response(200,
                ["Content-Type" => "application/json"],
                JSON3.write(summary))
        elseif fmt == "csv"
            return HTTP.Response(200,
                ["Content-Type" => "text/csv"],
                _export_csv(summary))
        else
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "format must be json or csv")))
        end
    end
end
