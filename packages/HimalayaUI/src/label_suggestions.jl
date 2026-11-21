using DBInterface, Tables

const ROLE_WORDS = (
    "Calibrant"       => Set(["agbe", "agbeh", "behenate", "calibrant", "lab6"]),
    "Empty capillary" => Set(["empty", "emptycap", "emptycapillary"]),
    "Buffer"          => Set(["buffer", "blank", "hepes", "pbs", "tris"]),
)

function _name_tokens(name::AbstractString)
    tokens = String[]
    for t in split(lowercase(name), r"[^a-z0-9]+"; keepempty = false)
        push!(tokens, t)
        bare = rstrip(isdigit, t)
        isempty(bare) || bare == t || push!(tokens, String(bare))
    end
    tokens
end

function suggested_role(name::AbstractString)
    tokens = _name_tokens(name)
    for (role, words) in ROLE_WORDS
        any(in(words), tokens) && return role
    end
    nothing
end

function suggest_labels!(db, sample_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT name FROM samples WHERE id = ? AND merged_into_id IS NULL", [sample_id]))
    (isempty(rows) || ismissing(rows[1].name)) && return nothing
    role = suggested_role(String(rows[1].name))
    role === nothing && return nothing
    roles = [first(r) for r in ROLE_WORDS]
    placeholders = join(fill("?", length(roles)), ", ")
    carried = Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM sample_labels WHERE sample_id = ? AND label COLLATE NOCASE IN ($placeholders) LIMIT 1",
        [sample_id; roles]))
    isempty(carried) || return nothing
    DBInterface.execute(db,
        """INSERT INTO sample_labels (sample_id, label, state, source) VALUES (?, ?, 'suggested', 'ingest')
           ON CONFLICT(sample_id, label COLLATE NOCASE) DO NOTHING""",
        [sample_id, role])
    role
end

function suggest_labels_for_experiment!(db, experiment_id::Int)
    ids = [Int(r.id) for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM samples WHERE experiment_id = ? AND merged_into_id IS NULL", [experiment_id]))]
    for id in ids
        suggest_labels!(db, id)
    end
    nothing
end
