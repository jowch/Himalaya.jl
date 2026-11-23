using DBInterface, Tables, SQLite

_same_dir(a, b) = _canonical_dir(a) == _canonical_dir(b)

function _setup_data_dir_taken(db::SQLite.DB, data_dir::AbstractString)
    plain = rstrip(String(data_dir), '/')
    !isempty(sqlrows(DBInterface.execute(db,
        "SELECT 1 FROM setups WHERE rtrim(data_dir, '/') IN (?, ?) LIMIT 1",
        [plain, _canonical_dir(data_dir)])))
end

function _unknown_setup_roots(db::SQLite.DB, experiment_id::Integer)::Vector{String}
    rows = sqlrows(DBInterface.execute(db,
        "SELECT path FROM experiments WHERE id = ?", [Int(experiment_id)]))
    isempty(rows) && return String[]
    path = String(first(rows).path)
    isdir(path) || return String[]
    known = setups_of(db, experiment_id)
    any(st -> _same_dir(st.root_dir, path), known) && return String[]
    roots = try
        setup_roots(path)
    catch err
        @warn "could not list setup folders" path exception = err
        return String[]
    end
    String[r for r in roots
           if !any(st -> _same_dir(st.root_dir, r), known) &&
              !_setup_data_dir_taken(db, _canonical_dir(_setup_data_dir(r)))]
end

_layout_detected(root) = resolve_experiment_layout(root).integration_pattern !== nothing

new_setup_roots(db::SQLite.DB, experiment_id::Integer) =
    filter(_layout_detected, _unknown_setup_roots(db, experiment_id))

unrecognized_setup_roots(db::SQLite.DB, experiment_id::Integer) =
    filter(!_layout_detected, _unknown_setup_roots(db, experiment_id))

function _next_position(db::SQLite.DB, experiment_id::Integer)
    rows = sqlrows(DBInterface.execute(db,
        "SELECT MAX(position) AS p FROM setups WHERE experiment_id = ?", [Int(experiment_id)]))
    p = first(rows).p
    ismissing(p) ? 0 : Int(p) + 1
end

function _folder_setup_spec(root::AbstractString)
    r = resolve_experiment_layout(root)
    data_dir = _canonical_dir(r.data_dir)
    (; data_dir,
       analysis_dir = something(r.analysis_dir, data_dir),
       root_dir = _canonical_dir(root),
       image_pattern = r.image_pattern,
       metadata_pattern = r.metadata_pattern,
       integration_pattern = r.integration_pattern)
end

function add_setup_from_folder!(db::SQLite.DB, experiment_id::Integer, root::AbstractString;
                                included::Bool = true)
    spec = _folder_setup_spec(root)
    create_setup!(db; spec..., experiment_id, included,
        position = _next_position(db, experiment_id))
end

function add_new_setups!(db::SQLite.DB, experiment_id::Integer)::Vector{Int}
    roots = new_setup_roots(db, experiment_id)
    isempty(roots) && return Int[]
    specs = [_folder_setup_spec(r) for r in roots]
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            ids = Int[]
            for spec in specs
                _setup_data_dir_taken(db, spec.data_dir) && continue
                push!(ids, create_setup!(db; spec..., experiment_id,
                    position = _next_position(db, experiment_id)))
            end
            ids
        end
    end
end

function _common_dir(paths::Vector{String})
    parts = [splitpath(rstrip(p, '/')) for p in paths]
    n = minimum(length, parts)
    k = 0
    while k < n && all(p -> p[k + 1] == parts[1][k + 1], parts)
        k += 1
    end
    k == 0 ? "/" : joinpath(parts[1][1:k]...)
end

function _experiment_summary(db::SQLite.DB, id::Integer)
    rows = sqlrows(DBInterface.execute(db, """
        SELECT x.id, x.name, x.path, x.ingest_status,
               (SELECT COUNT(*) FROM samples s WHERE s.experiment_id = x.id) AS samples,
               (SELECT COUNT(*) FROM exposures e WHERE e.experiment_id = x.id) AS exposures
          FROM experiments x WHERE x.id = ?""", [Int(id)]))
    isempty(rows) && error("No experiment with id $id.")
    first(rows)
end

function merge_plan(db::SQLite.DB, ids::Vector{Int};
                    name::Union{Nothing,AbstractString} = nothing,
                    path::Union{Nothing,AbstractString} = nothing)
    length(unique(ids)) == length(ids) || error("An experiment is listed twice.")
    length(ids) >= 2 || error("Name at least two experiments to merge.")
    exps = [_experiment_summary(db, id) for id in ids]
    for x in exps
        x.ingest_status in ("scanning", "analyzing") &&
            error("Experiment $(x.id) ($(x.name)) is being scanned. Wait for it to finish.")
    end
    setups = reduce(vcat, [collect(setups_of(db, id)) for id in ids])
    isempty(setups) && error("None of these experiments has a setup.")
    folders = [st.folder for st in setups]
    dupes = unique(f for f in folders if count(==(f), folders) > 1)
    isempty(dupes) || error("Two setups share the folder name $(join(dupes, ", ")). " *
                            "Merge only setups from one beamtime folder.")
    roots = String[String(st.root_dir) for st in setups]
    if path === nothing
        parents = unique(_canonical_dir(dirname(rstrip(r, '/'))) for r in roots)
        length(parents) == 1 || error("These setups are not in one beamtime folder. " *
                                      "Pass --path to choose the beamtime folder.")
    end
    merged_path = path === nothing ? _common_dir(roots) : _canonical_dir(String(path))
    holder = experiment_at_path(db, merged_path)
    holder === nothing || holder.id in ids ||
        error("Experiment $(holder.id) ($(holder.name)) already uses $(merged_path).")
    merged_name = name === nothing ? basename(merged_path) : String(name)
    order = sortperm(folders)
    (; into = first(ids), from = ids[2:end], experiments = exps,
       name = merged_name, path = merged_path,
       setups = [(; id = Int(setups[i].id), folder = String(setups[i].folder),
                    experiment_id = Int(setups[i].experiment_id), position = p - 1)
                 for (p, i) in enumerate(order)])
end

function apply_merge!(db::SQLite.DB, plan)
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            for st in plan.setups
                DBInterface.execute(db, "UPDATE setups SET experiment_id = ?, position = ? WHERE id = ?",
                    [plan.into, st.position, st.id])
            end
            for from in plan.from, t in ("loads", "samples", "exposures")
                DBInterface.execute(db, "UPDATE $t SET experiment_id = ? WHERE experiment_id = ?",
                    [plan.into, from])
            end
            for from in plan.from
                DBInterface.execute(db, "DELETE FROM experiments WHERE id = ?", [from])
            end
            DBInterface.execute(db, "UPDATE experiments SET name = ?, path = ? WHERE id = ?",
                [plan.name, plan.path, plan.into])
        end
    end
    nothing
end

function print_merge_plan(io::IO, plan)
    println(io, "Merge into experiment $(plan.into):")
    for x in plan.experiments
        println(io, "  [$(x.id)] $(x.name)  ($(x.path))  $(x.samples) samples, $(x.exposures) images")
    end
    println(io, "Result: \"$(plan.name)\" at $(plan.path)")
    for st in plan.setups
        println(io, "  setup $(st.position + 1): $(st.folder)  (from experiment $(st.experiment_id))")
    end
    nothing
end
