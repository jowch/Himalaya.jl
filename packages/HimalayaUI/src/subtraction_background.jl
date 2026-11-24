using DBInterface, SQLite

const BUFFER_SUFFIX = r"_B\d+$"

is_buffer_filename(filename::AbstractString) = occursin(AUTOSAMPLER_STEM, filename) && occursin(BUFFER_SUFFIX, filename)

function autosampler_runs(db::SQLite.DB, experiment_id::Integer)
    rows = sqlrows(DBInterface.execute(db, """
        SELECT s.id AS id, s.load_id AS load_id, l.load_index AS load_index, s.slot_index AS slot, e.filename AS filename
          FROM samples s
          JOIN setups st ON st.id = s.setup_id
          LEFT JOIN loads l ON l.id = s.load_id
          LEFT JOIN exposures e ON e.sample_id = s.id
         WHERE s.experiment_id = ? AND st.mode = 'autosampler'
         ORDER BY s.id, e.id""", [Int(experiment_id)]))
    load_of = Dict{Int,Union{Int,Missing}}()
    slot_of = Dict{Int,Int}()
    buffers = Set{Int}()
    for r in rows
        id = Int(r.id)
        in_run = !ismissing(r.load_id) && !ismissing(r.load_index) && r.load_index != 0
        load_of[id] = in_run ? Int(r.load_id) : missing
        slot_of[id] = ismissing(r.slot) ? typemax(Int) : Int(r.slot)
        !ismissing(r.filename) && is_buffer_filename(String(r.filename)) && push!(buffers, id)
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

default_background_of(run) = run.is_buffer || isempty(run.run_buffers) ? nothing : first(run.run_buffers)

default_backgrounds(runs::AbstractDict) =
    Dict{Int,Union{Int,Nothing}}(id => default_background_of(run) for (id, run) in runs)

default_backgrounds(db::SQLite.DB, experiment_id::Integer) = default_backgrounds(autosampler_runs(db, experiment_id))

function default_background(db::SQLite.DB, sample_id::Integer)
    rows = sqlrows(DBInterface.execute(db, "SELECT experiment_id FROM samples WHERE id = ?", [Int(sample_id)]))
    isempty(rows) || ismissing(rows[1].experiment_id) ? nothing :
        get(default_backgrounds(db, Int(rows[1].experiment_id)), Int(sample_id), nothing)
end
