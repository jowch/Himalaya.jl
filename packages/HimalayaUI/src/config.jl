using TOML

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
        get(mf,  "sample_id",      1),
        get(mf,  "label",          2),
        get(mf,  "name",           3),
        get(mf,  "filenames",      9),
        get(mf,  "notes_sample",   10),
        get(mf,  "notes_exposure", 11),
        get(layout, "data_dir",      "data"),
        get(layout, "analysis_dir",  "analysis/automatic_analysis"),
        get(layout, "exposure_type", "simple"),
        integration,
        image,
    )
end
