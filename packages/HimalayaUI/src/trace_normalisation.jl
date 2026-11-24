using DBInterface, SQLite

const INTEG_FILE = "integ.mpp"
const NORMALISED_BY_I2 = "i2"
const NORMALISED_BY_I1 = "i1"
const RAW_COUNTS = "none"
const NORMALISATION_UNKNOWN = "unknown"

_dir_key(d) = rstrip(normpath(abspath(String(d))), '/')
_inside(dir, root) = (d = _dir_key(dir); r = _dir_key(root); d == r || startswith(d, r * "/"))

function find_integ_file(analysis_dir::AbstractString; within::AbstractString = analysis_dir)::Union{String,Nothing}
    dir = _dir_key(analysis_dir)
    while true
        candidate = joinpath(dir, INTEG_FILE)
        isfile(candidate) && return candidate
        parent = dirname(dir)
        (parent == dir || !_inside(parent, within)) && return nothing
        dir = parent
    end
end

# physics: SASTool's `-i yes <i1|i2> <grand base>` divides every trace by that
# physics: monitor, times a grand base shared by the series. Without the line the
# physics: traces are raw counts. Which one applies decides the buffer scale.
function parse_integ_normalisation(path::AbstractString)
    for line in eachline(path)
        tokens = split(strip(line))
        (isempty(tokens) || tokens[1] != "-i") && continue
        length(tokens) >= 2 && lowercase(tokens[2]) == "yes" || return (by = RAW_COUNTS, base = nothing)
        length(tokens) >= 3 || return (by = NORMALISATION_UNKNOWN, base = nothing)
        channel = lowercase(tokens[3])
        channel in (NORMALISED_BY_I1, NORMALISED_BY_I2) || return (by = NORMALISATION_UNKNOWN, base = nothing)
        base = length(tokens) >= 4 ? tryparse(Float64, tokens[4]) : nothing
        return (by = channel, base = base)
    end
    (by = RAW_COUNTS, base = nothing)
end

function trace_normalisation(analysis_dir::AbstractString; within::AbstractString = analysis_dir)
    path = find_integ_file(analysis_dir; within)
    path === nothing && return (by = NORMALISATION_UNKNOWN, base = nothing)
    try
        parse_integ_normalisation(path)
    catch err
        @warn "could not read the integration settings" path exception = err
        (by = NORMALISATION_UNKNOWN, base = nothing)
    end
end

function read_trace_normalisations!(db::SQLite.DB; experiment_id::Union{Integer,Nothing} = nothing)::Int
    sql = "SELECT s.id, s.analysis_dir, s.root_dir, e.path, s.trace_normalisation, s.normalisation_base
             FROM setups s JOIN experiments e ON e.id = s.experiment_id" *
          (experiment_id === nothing ? "" : " WHERE s.experiment_id = ?")
    rows = sqlrows(DBInterface.execute(db, sql, experiment_id === nothing ? [] : [Int(experiment_id)]))
    changed = Tuple{Int,NamedTuple}[]
    for r in rows
        isdir(String(r.analysis_dir)) || continue
        root = ismissing(r.path) || isempty(String(r.path)) ? String(r.root_dir) : String(r.path)
        n = trace_normalisation(String(r.analysis_dir); within = root)
        isequal(coalesce(r.trace_normalisation, nothing), n.by) &&
            isequal(coalesce(r.normalisation_base, nothing), n.base) && continue
        push!(changed, (Int(r.id), n))
    end
    isempty(changed) && return 0
    lock(_DB_WRITE_LOCK) do
        for (id, n) in changed
            DBInterface.execute(db,
                "UPDATE setups SET trace_normalisation = ?, normalisation_base = ? WHERE id = ?",
                [n.by, n.base, id])
        end
    end
    length(changed)
end
