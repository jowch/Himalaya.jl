using DBInterface, SQLite

const SETUP_MODES = ("capillary", "autosampler")
const MODE_SOURCES = ("default", "detected", "user", "migration")

const DEFAULT_FRAME_SUFFIX = "_0_001"
const AUTOSAMPLER_INTEGRATION_PATTERN = "{name}_tot.dat"

const AUTOSAMPLER_RUN_FILE = r"^AT_SolSAXS_.+_(\d+)\.txt$"
const WELL_TOKEN = r"_\d{2}[A-Z]_[SB]\d+(_\d+_\d+)?$"

is_autosampler(setup) =
    hasproperty(setup, :mode) && !ismissing(setup.mode) && setup.mode == "autosampler"

function detect_setup_mode(data_dir::AbstractString)::Union{String,Nothing}
    isdir(data_dir) || return nothing
    names = list_names(data_dir)
    any(f -> occursin(AUTOSAMPLER_RUN_FILE, f), names) && return "autosampler"
    any(f -> endswith(f, ".tif") && !startswith(f, "."), names) ? "capillary" : nothing
end

function first_frame_suffix(data_dir::AbstractString)::String
    isdir(data_dir) || return DEFAULT_FRAME_SUFFIX
    names = list_names(data_dir)
    best = nothing
    for f in names
        startswith(f, ".") && continue
        m = match(r"(_\d+_(\d+))\.tif$", f)
        m === nothing && continue
        key = (parse(Int, m.captures[2]), String(m.captures[1]))
        (best === nothing || key < best) && (best = key)
    end
    best === nothing ? DEFAULT_FRAME_SUFFIX : best[2]
end

const AUTOSAMPLER_NOT_EDITABLE = "An autosampler sample is not analyzed, so its peaks, calls and indices can't be edited."
const AUTOSAMPLER_NOT_SELECTABLE = "An autosampler sample has one image, so there is nothing to select."

_autosampler_refusal(message) = HTTP.Response(422, ["Content-Type" => "application/json"],
    to_json(Dict(:error => message)))

_in_autosampler_setup(db::SQLite.DB, sql::AbstractString, id::Integer) =
    !isempty(sqlrows(DBInterface.execute(db, sql, [Int(id)])))

function refuse_autosampler_sample(db::SQLite.DB, sample_id::Integer, message::AbstractString = AUTOSAMPLER_NOT_EDITABLE)
    _in_autosampler_setup(db, """
        SELECT 1 AS x FROM samples s JOIN setups st ON st.id = s.setup_id
         WHERE s.id = ? AND st.mode = 'autosampler'""", sample_id) ? _autosampler_refusal(message) : nothing
end

function refuse_autosampler_exposure(db::SQLite.DB, exposure_id::Integer)
    _in_autosampler_setup(db, """
        SELECT 1 AS x FROM exposures e JOIN samples s ON s.id = e.sample_id
          JOIN setups st ON st.id = s.setup_id
         WHERE e.id = ? AND st.mode = 'autosampler'""", exposure_id) ? _autosampler_refusal(AUTOSAMPLER_NOT_EDITABLE) : nothing
end

const AUTOSAMPLER_NOT_BRIDGEABLE = "A capillary sample and an autosampler sample can't be the same sample."

function autosampler_samples(db::SQLite.DB, sample_ids::AbstractVector{<:Integer})
    isempty(sample_ids) && return Set{Int}()
    marks = join(fill("?", length(sample_ids)), ", ")
    Set(Int(r.id) for r in sqlrows(DBInterface.execute(db, """
        SELECT s.id FROM samples s JOIN setups st ON st.id = s.setup_id
         WHERE st.mode = 'autosampler' AND s.id IN ($marks)""", Int.(collect(sample_ids)))))
end

function refuse_autosampler_row(db::SQLite.DB, table::AbstractString, id::Integer)
    rows = sqlrows(DBInterface.execute(db, "SELECT exposure_id FROM $table WHERE id = ?", [Int(id)]))
    isempty(rows) || ismissing(rows[1].exposure_id) ? nothing : refuse_autosampler_exposure(db, Int(rows[1].exposure_id))
end

function autosampler_layout(root::AbstractString)
    suffix = first_frame_suffix(_setup_data_dir(String(root)))
    (; dir = joinpath(String(root), "analysis", "automatic_analysis", "tot_files"),
       image_pattern = "{name}$(suffix).tif",
       metadata_pattern = "{name}$(suffix).prp",
       integration_pattern = AUTOSAMPLER_INTEGRATION_PATTERN)
end

_setup_has_exposures(db::SQLite.DB, setup_id::Int) =
    !isempty(sqlrows(DBInterface.execute(db,
        "SELECT 1 FROM exposures WHERE setup_id = ? LIMIT 1", [setup_id])))

_has_well_token_exposures(db::SQLite.DB, setup_id::Int) =
    any(r -> occursin(WELL_TOKEN, String(r.filename)), sqlrows(DBInterface.execute(db,
        "SELECT filename FROM exposures WHERE setup_id = ? AND filename IS NOT NULL", [setup_id])))

_has_reviewed_exposures(db::SQLite.DB, setup_id::Int) =
    !isempty(sqlrows(DBInterface.execute(db,
        "SELECT 1 FROM exposures WHERE setup_id = ? AND reviewed_at IS NOT NULL LIMIT 1", [setup_id])))

function _may_promote(db::SQLite.DB, setup)::Bool
    open_to_detection = ismissing(setup.mode) ||
        (setup.mode == "capillary" && coalesce(setup.mode_source, "default") == "detected")
    open_to_detection && !_has_well_token_exposures(db, Int(setup.id)) &&
        !_has_reviewed_exposures(db, Int(setup.id))
end

function _write_layout!(db::SQLite.DB, setup_id::Int, layout)
    DBInterface.execute(db, """
        UPDATE setups SET analysis_dir = ?, image_pattern = ?, metadata_pattern = ?,
                          integration_pattern = ? WHERE id = ?""",
        [layout.dir, layout.image_pattern, layout.metadata_pattern, layout.integration_pattern, setup_id])
end

function _promote_to_autosampler!(db::SQLite.DB, setup)
    id = Int(setup.id)
    layout = autosampler_layout(String(setup.root_dir))
    first_frame = replace(layout.image_pattern, "{name}" => "", ".tif" => "")
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            DBInterface.execute(db,
                "UPDATE setups SET mode = 'autosampler', mode_source = 'detected' WHERE id = ?", [id])
            _write_layout!(db, id, layout)
            for r in sqlrows(DBInterface.execute(db,
                    "SELECT id, filename FROM exposures WHERE setup_id = ? AND filename IS NOT NULL", [id]))
                name = String(r.filename)
                endswith(name, first_frame) && DBInterface.execute(db,
                    "UPDATE exposures SET filename = ? WHERE id = ?",
                    [chop(name; tail = length(first_frame)), Int(r.id)])
            end
        end
    end
    setup_row(db, id)
end

function _relayout_empty_autosampler!(db::SQLite.DB, setup)
    id = Int(setup.id)
    _setup_has_exposures(db, id) && return setup
    layout = autosampler_layout(String(setup.root_dir))
    same = !ismissing(setup.analysis_dir) && !ismissing(setup.image_pattern) &&
           !ismissing(setup.metadata_pattern) && !ismissing(setup.integration_pattern) &&
           (setup.analysis_dir, setup.image_pattern, setup.metadata_pattern, setup.integration_pattern) ==
           (layout.dir, layout.image_pattern, layout.metadata_pattern, layout.integration_pattern)
    same && return setup
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            _write_layout!(db, id, layout)
        end
    end
    setup_row(db, id)
end

function settle_setup_mode!(db::SQLite.DB, setup)
    is_autosampler(setup) && return _relayout_empty_autosampler!(db, setup)
    detected = detect_setup_mode(String(setup.data_dir))
    detected === nothing && return setup
    detected == "autosampler" && _may_promote(db, setup) && return _promote_to_autosampler!(db, setup)
    ismissing(setup.mode) || return setup
    lock(_DB_WRITE_LOCK) do
        DBInterface.execute(db,
            "UPDATE setups SET mode = 'capillary', mode_source = 'detected' WHERE id = ?", [Int(setup.id)])
    end
    setup_row(db, Int(setup.id))
end

function is_autosampler_sample(db::SQLite.DB, sample_id::Integer)::Bool
    rows = sqlrows(DBInterface.execute(db, """
        SELECT st.mode FROM samples s JOIN setups st ON st.id = s.setup_id WHERE s.id = ?""", [Int(sample_id)]))
    !isempty(rows) && is_autosampler(first(rows))
end
