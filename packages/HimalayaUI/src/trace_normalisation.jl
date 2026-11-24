using DBInterface, SQLite

const INTEG_FILE = "integ.mpp"
const NORMALISED_BY_I2 = "i2"
const NORMALISED_BY_I1 = "i1"
const RAW_COUNTS = "none"
const NORMALISATION_UNKNOWN = "unknown"

function find_integ_file(analysis_dir::AbstractString)::Union{String,Nothing}
    dir = String(analysis_dir)
    for _ in 1:4
        candidate = joinpath(dir, INTEG_FILE)
        isfile(candidate) && return candidate
        parent = dirname(dir)
        parent == dir && break
        dir = parent
    end
    nothing
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

function trace_normalisation(analysis_dir::AbstractString)
    path = find_integ_file(analysis_dir)
    path === nothing && return (by = NORMALISATION_UNKNOWN, base = nothing)
    try
        parse_integ_normalisation(path)
    catch err
        @warn "could not read the integration settings" path exception = err
        (by = NORMALISATION_UNKNOWN, base = nothing)
    end
end

function read_trace_normalisations!(db::SQLite.DB; experiment_id::Union{Integer,Nothing} = nothing)::Int
    sql = "SELECT id, analysis_dir FROM setups
            WHERE (trace_normalisation IS NULL OR trace_normalisation = '$NORMALISATION_UNKNOWN')" *
          (experiment_id === nothing ? "" : " AND experiment_id = ?")
    rows = sqlrows(DBInterface.execute(db, sql, experiment_id === nothing ? [] : [Int(experiment_id)]))
    found = [(Int(r.id), trace_normalisation(String(r.analysis_dir))) for r in rows]
    isempty(found) && return 0
    lock(_DB_WRITE_LOCK) do
        for (id, n) in found
            DBInterface.execute(db,
                "UPDATE setups SET trace_normalisation = ?, normalisation_base = ? WHERE id = ?",
                [n.by, n.base, id])
        end
    end
    count(f -> f[2].by != NORMALISATION_UNKNOWN, found)
end
