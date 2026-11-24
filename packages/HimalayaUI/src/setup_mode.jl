using DBInterface, SQLite

const SETUP_MODES = ("capillary", "autosampler")

const AUTOSAMPLER_IMAGE_PATTERN = "{name}_0_001.tif"
const AUTOSAMPLER_METADATA_PATTERN = "{name}_0_001.prp"
const AUTOSAMPLER_INTEGRATION_PATTERN = "{name}_tot.dat"

const WELL_TOKEN = r"_\d{2}[A-Z]_[SB]\d+(_0_\d+)?$"

is_autosampler(setup) =
    hasproperty(setup, :mode) && !ismissing(setup.mode) && setup.mode == "autosampler"

function detect_setup_mode(data_dir::AbstractString)::Union{String,Nothing}
    names = try
        readdir(String(data_dir))
    catch
        return nothing
    end
    any(f -> occursin(AUTOSAMPLER_RUN_FILE, f), names) && return "autosampler"
    any(f -> endswith(f, ".tif") && !startswith(f, "."), names) ? "capillary" : nothing
end

function autosampler_layout(root::AbstractString)
    tot = joinpath(String(root), "analysis", "automatic_analysis", "tot_files")
    (; dir = isdir(tot) ? tot : nothing,
       image_pattern = AUTOSAMPLER_IMAGE_PATTERN,
       metadata_pattern = AUTOSAMPLER_METADATA_PATTERN,
       integration_pattern = AUTOSAMPLER_INTEGRATION_PATTERN)
end

function settle_setup_mode!(db::SQLite.DB, setup)
    is_autosampler(setup) && return setup
    detected = detect_setup_mode(String(setup.data_dir))
    detected === nothing && return setup
    promote = detected == "autosampler" && !_has_well_token_exposures(db, Int(setup.id))
    new_mode = promote ? "autosampler" : (ismissing(setup.mode) ? "capillary" : nothing)
    new_mode === nothing && return setup
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            DBInterface.execute(db, "UPDATE setups SET mode = ? WHERE id = ?", [new_mode, Int(setup.id)])
            if promote
                layout = autosampler_layout(String(setup.root_dir))
                DBInterface.execute(db, """
                    UPDATE setups SET analysis_dir = ?, image_pattern = ?, metadata_pattern = ?,
                                      integration_pattern = ? WHERE id = ?""",
                    [something(layout.dir, String(setup.analysis_dir)), layout.image_pattern,
                     layout.metadata_pattern, layout.integration_pattern, Int(setup.id)])
            end
        end
    end
    setup_row(db, Int(setup.id))
end

_has_well_token_exposures(db::SQLite.DB, setup_id::Int) =
    any(r -> occursin(WELL_TOKEN, String(r.filename)), sqlrows(DBInterface.execute(db,
        "SELECT filename FROM exposures WHERE setup_id = ? AND filename IS NOT NULL", [setup_id])))

function is_autosampler_sample(db::SQLite.DB, sample_id::Integer)::Bool
    rows = sqlrows(DBInterface.execute(db, """
        SELECT st.mode FROM samples s JOIN setups st ON st.id = s.setup_id WHERE s.id = ?""", [Int(sample_id)]))
    !isempty(rows) && is_autosampler(first(rows))
end
