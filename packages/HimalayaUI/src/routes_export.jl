using HTTP, JSON3, DBInterface, Tables, Oxygen, CSV, SQLite

"""
    build_export(db, experiment_id) -> Vector{Dict}

Per-sample summary with the sample's call; `call_state` is `uncalled` until someone ticks a candidate.
A called sample lists each exposure's members (none when that frame adds nothing, still `called`); others list every candidate.
"""
function build_export(db::SQLite.DB, experiment_id::Int)
    samples = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM samples WHERE experiment_id = ? ORDER BY id",
        [experiment_id]))
    sids = [Int(sm.id) for sm in samples]
    states = sample_call_states(db, sids)
    members = sample_call_members(db, sids)

    out = Dict[]
    for sm in samples
        sid  = Int(sm.id)
        labels = [String(r.label) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT label FROM sample_labels WHERE sample_id = ? AND state = 'confirmed' ORDER BY id", [sid]))]
        state = get(states, sid, "indexed")
        ms = get(members, sid, [])
        summary = call_summary(state, ms)
        called = summary.phase !== nothing
        called_ids = Set(m.index_id for m in ms)

        exposures = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM exposures WHERE sample_id = ? ORDER BY id", [sid]))

        ex_out = Dict[]
        for e in exposures
            indices = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, phase, basis, score, r_squared, lattice_d
                 FROM indices WHERE exposure_id = ? ORDER BY score DESC, id", [Int(e.id)]))
            called && (indices = filter(ix -> Int(ix.id) in called_ids, indices))
            push!(ex_out, Dict(
                :filename => e.filename,
                :selected => e.selected != 0,
                :called   => called,
                :indices  => [delete!(ix, :id) for ix in rows_to_json(indices)],
            ))
        end

        push!(out, Dict(
            :name        => sm.name,
            :notes       => sm.notes isa Missing ? nothing : sm.notes,
            :labels      => labels,
            :call_state  => export_call_state(state, ms),
            :phase       => summary.phase,
            :call_phases => summary.call_phases,
            :exposures   => ex_out,
        ))
    end
    out
end

function _export_csv(summary)
    header = ["sample_name","exposure_filename","call_state","call_phases",
              "phases","lattice_ds","r_squareds","scores","labels"]
    rows = Vector{Vector{Any}}()
    for sm in summary
        label_str = join(sm[:labels], ";")
        call_str = join(sm[:call_phases], ";")
        if isempty(sm[:exposures])
            push!(rows, [sm[:name], "", sm[:call_state], call_str, "", "", "", "", label_str])
            continue
        end
        for ex in sm[:exposures]
            phases = join([String(ix[:phase])     for ix in ex[:indices]], ";")
            lds    = join([string(ix[:lattice_d]) for ix in ex[:indices]], ";")
            r2s    = join([string(ix[:r_squared]) for ix in ex[:indices]], ";")
            scs    = join([string(ix[:score])     for ix in ex[:indices]], ";")
            push!(rows, [sm[:name], ex[:filename], sm[:call_state], call_str,
                         phases, lds, r2s, scs, label_str])
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
