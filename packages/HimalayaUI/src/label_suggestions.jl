using DBInterface, Tables, SQLite

const BUFFER_WORDS = Set(["buffer", "buf", "blank"])

const ROLE_WORDS = (
    "Calibrant"       => Set(["agbe", "agbeh", "behenate", "calibrant", "lab6"]),
    "Empty capillary" => Set(["empty", "emptycap", "emptycapillary"]),
    "Buffer"          => BUFFER_WORDS,
)

const SOLVENT_WORDS = Set(["hepes", "pbs", "tris"])

_bare(t) = (b = rstrip(isdigit, t); isempty(b) ? t : String(b))
_is_code(t) = occursin(r"^[a-z]?[0-9]+$", t)

function suggested_role(name::AbstractString)
    tokens = split(lowercase(name), r"[^a-z0-9]+"; keepempty = false)
    words = [w for t in tokens for w in unique((String(t), _bare(t)))]
    for (role, vocab) in ROLE_WORDS
        any(in(vocab), words) && return role
    end
    any(t -> _bare(t) in SOLVENT_WORDS, tokens) &&
        all(t -> _is_code(t) || _bare(t) in SOLVENT_WORDS || _bare(t) in BUFFER_WORDS, tokens) &&
        return "Buffer"
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

function suggest_labels_quietly!(db, experiment_id::Int)
    try
        lock(_DB_WRITE_LOCK) do
            SQLite.transaction(() -> suggest_labels_for_experiment!(db, experiment_id), db)
        end
    catch err
        @warn "label suggestions failed" experiment_id exception = err
    end
    nothing
end

function migrate_label_suggestions!(db::SQLite.DB)
    _migrated(db, MIGRATION_LABEL_SUGGESTIONS) && return nothing
    SQLite.transaction(db) do
        for r in Tables.rowtable(DBInterface.execute(db, "SELECT id FROM experiments"))
            suggest_labels_for_experiment!(db, Int(r.id))
        end
        _record_migration!(db, MIGRATION_LABEL_SUGGESTIONS)
    end
    nothing
end
