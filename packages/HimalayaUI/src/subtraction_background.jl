using DBInterface, SQLite

const BUFFER_FILENAME = r"_B\d+$"

function autosampler_runs(db::SQLite.DB, experiment_id::Integer)
    rows = sqlrows(DBInterface.execute(db, """
        SELECT s.id AS id, s.load_id AS load_id, s.slot_index AS slot, e.filename AS filename
          FROM samples s
          JOIN setups st ON st.id = s.setup_id
          LEFT JOIN exposures e ON e.sample_id = s.id
         WHERE s.experiment_id = ? AND st.mode = 'autosampler'
         ORDER BY s.id, e.id""", [Int(experiment_id)]))
    load_of = Dict{Int,Union{Int,Missing}}()
    slot_of = Dict{Int,Int}()
    buffers = Set{Int}()
    for r in rows
        id = Int(r.id)
        load_of[id] = ismissing(r.load_id) ? missing : Int(r.load_id)
        slot_of[id] = ismissing(r.slot) ? typemax(Int) : Int(r.slot)
        !ismissing(r.filename) && occursin(BUFFER_FILENAME, String(r.filename)) && push!(buffers, id)
    end
    run_buffers = Dict{Int,Vector{Int}}()
    for id in sort(collect(buffers); by = id -> (slot_of[id], id))
        load = load_of[id]
        ismissing(load) || push!(get!(run_buffers, load, Int[]), id)
    end
    Dict(id => (; is_buffer = id in buffers,
                  run_buffers = ismissing(load) ? Int[] : get(run_buffers, load, Int[]))
         for (id, load) in load_of)
end

function default_backgrounds(db::SQLite.DB, experiment_id::Integer)
    Dict{Int,Union{Int,Nothing}}(
        id => (run.is_buffer || isempty(run.run_buffers) ? nothing : first(run.run_buffers))
        for (id, run) in autosampler_runs(db, experiment_id))
end

function default_background(db::SQLite.DB, sample_id::Integer)
    rows = sqlrows(DBInterface.execute(db, "SELECT experiment_id FROM samples WHERE id = ?", [Int(sample_id)]))
    isempty(rows) || ismissing(rows[1].experiment_id) ? nothing :
        get(default_backgrounds(db, Int(rows[1].experiment_id)), Int(sample_id), nothing)
end
