using TOML

"""
    ExperimentConfig

Parsed representation of an `experiment.toml` file. Drives all file I/O conventions
for an experiment: manifest column layout, directory paths, file name patterns,
beamline parameters.

Supported TOML sections:
- `[experiment]`: name, description, manifest (relative path to manifest CSV)
- `[beamline]`: energy_kev, flight_path_m
- `[manifest]`: delimiter, skip_rows, header_row, sample_id, label, name,
  filenames, notes_sample, notes_exposure (each column = Int index or String header name)
- `[layout]`: data_dir, analysis_dir, exposure_type
- `[files]`: integration, image (path patterns containing `{name}`)
"""
struct ExperimentConfig
    # [experiment]
    name               ::String
    description        ::String
    manifest_file      ::String
    # [beamline]
    energy_kev         ::Union{Float64,Nothing}
    flight_path_m      ::Union{Float64,Nothing}
    # [manifest]
    delimiter          ::String
    skip_rows          ::Int
    header_row         ::Int
    col_sample_id      ::Union{Int,String}
    col_label          ::Union{Int,String}
    col_name           ::Union{Int,String}
    col_filenames      ::Union{Int,String}
    col_notes_sample   ::Union{Int,String}
    col_notes_exposure ::Union{Int,String}
    # [layout]
    data_dir           ::String
    analysis_dir       ::String
    exposure_type      ::String
    # [files]
    integration_pattern::String
    image_pattern      ::String
end

function _validate_pattern(pattern::String, field::String)
    startswith(pattern, "/") && error("$field must not be an absolute path: $pattern")
    contains(pattern, "..") && error("$field must not traverse upward (..): $pattern")
    contains(pattern, "{name}") || error("$field must contain {name}: $pattern")
end

function _coerce_col(v, field::String)::Union{Int,String}
    v isa Integer        && return Int(v)
    v isa AbstractString && return String(v)
    error("$field must be an integer column index or a string header name, got $(typeof(v)): $v")
end

"""
    load_config(path) -> ExperimentConfig

Read and parse an `experiment.toml` file at `path`. Validates that file-name
patterns in `[files]` are relative, do not traverse upward (`..`), and contain
`{name}`. Validates that each `[manifest]` column entry is either an integer
index or a string header name. Throws `ErrorException` on any violation, or if
the file does not exist.
"""
function load_config(path::AbstractString)::ExperimentConfig
    isfile(path) || error("experiment.toml not found: $path")
    d = TOML.parsefile(path)

    exp    = get(d, "experiment", Dict())
    bl     = get(d, "beamline",   Dict())
    mf     = get(d, "manifest",   Dict())
    layout = get(d, "layout",     Dict())
    files  = get(d, "files",      Dict())

    integration = get(files, "integration", "{name}.dat")
    image       = get(files, "image",       "{name}.tiff")
    _validate_pattern(integration, "files.integration")
    _validate_pattern(image,       "files.image")

    ExperimentConfig(
        get(exp, "name",        ""),
        get(exp, "description", ""),
        get(exp, "manifest",    "manifest.csv"),
        get(bl,  "energy_kev",    nothing),
        get(bl,  "flight_path_m", nothing),
        get(mf,  "delimiter",      "\t"),
        get(mf,  "skip_rows",      1),
        get(mf,  "header_row",     0),
        _coerce_col(get(mf, "sample_id",      1),  "manifest.sample_id"),
        _coerce_col(get(mf, "label",          2),  "manifest.label"),
        _coerce_col(get(mf, "name",           3),  "manifest.name"),
        _coerce_col(get(mf, "filenames",      9),  "manifest.filenames"),
        _coerce_col(get(mf, "notes_sample",   10), "manifest.notes_sample"),
        _coerce_col(get(mf, "notes_exposure", 11), "manifest.notes_exposure"),
        get(layout, "data_dir",      "data"),
        get(layout, "analysis_dir",  "analysis/automatic_analysis"),
        get(layout, "exposure_type", "simple"),
        integration,
        image,
    )
end
