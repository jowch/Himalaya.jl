using TOML

const VALID_EXPOSURE_TYPES = ("simple",)

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
    d = try
        TOML.parsefile(path)
    catch e
        error("Invalid TOML in $path: $(sprint(showerror, e))")
    end
    _build_config(d)
end

function _build_config(d::AbstractDict)::ExperimentConfig
    exp    = get(d, "experiment", Dict())
    bl     = get(d, "beamline",   Dict())
    mf     = get(d, "manifest",   Dict())
    layout = get(d, "layout",     Dict())
    files  = get(d, "files",      Dict())

    integration = get(files, "integration", "{name}.dat")
    image       = get(files, "image",       "{name}.tiff")
    _validate_pattern(integration, "files.integration")
    _validate_pattern(image,       "files.image")

    exposure_type = get(layout, "exposure_type", "simple")
    exposure_type in VALID_EXPOSURE_TYPES || error(
        "layout.exposure_type '$exposure_type' not recognized. Valid: $(join(VALID_EXPOSURE_TYPES, ", "))")

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
        exposure_type,
        integration,
        image,
    )
end

"""
    configs_dir() -> String

Returns the path to the directory containing config templates. When
`HIMALAYA_CONFIGS_DIR` is set in the environment, that path is used
(allowing labs to layer in custom templates without patching the
package). Otherwise returns the package's bundled `configs/` directory.
"""
configs_dir() = get(ENV, "HIMALAYA_CONFIGS_DIR", joinpath(@__DIR__, "..", "configs"))

"""
    list_config_types() -> Vector{String}

Returns the names of all built-in config templates (filename stems of
`*.toml` files in `configs_dir()`).
"""
function list_config_types()::Vector{String}
    dir = configs_dir()
    isdir(dir) || return String[]
    [splitext(f)[1] for f in readdir(dir) if endswith(f, ".toml")]
end

"""
    load_builtin_config(type_name) -> ExperimentConfig

Loads a built-in config template by name (without `.toml` extension).
Errors with the list of available types if the name is unknown.
"""
function load_builtin_config(type_name::AbstractString)::ExperimentConfig
    path = joinpath(configs_dir(), type_name * ".toml")
    isfile(path) || error("Unknown config type '$type_name'. Available: $(join(list_config_types(), ", "))")
    load_config(path)
end

"""
    resolve_files(cfg, base_dir, prefix, file_pattern) -> Vector{String}

Scan `base_dir` for files matching `file_pattern` where `{name}` is replaced by
`prefix*` (any string starting with `prefix`). Returns sorted bare filename stems
with the pattern's trailing suffix stripped (e.g. "JC001" not "JC001.dat").

`file_pattern` may contain a leading subdirectory (e.g. `"integrated/{name}.dat"`),
in which case scanning happens in `joinpath(base_dir, "integrated")`. Returns
an empty vector if the scan directory doesn't exist.
"""
function resolve_files(
    ::ExperimentConfig,
    base_dir::AbstractString,
    prefix::AbstractString,
    file_pattern::String,
)::Vector{String}
    parts = split(file_pattern, "{name}"; limit=2)
    length(parts) == 2 || error("file pattern must contain exactly one {name}: $file_pattern")
    before, after = String(parts[1]), String(parts[2])

    scan_subdir       = dirname(before)
    file_prefix_local = basename(before) * prefix
    scan_dir = isempty(scan_subdir) ? String(base_dir) : joinpath(base_dir, scan_subdir)

    if !isdir(scan_dir)
        @warn "resolve_files: scan directory does not exist" scan_dir pattern=file_pattern
        return String[]
    end

    matches = filter(readdir(scan_dir)) do f
        startswith(f, file_prefix_local) && endswith(f, after)
    end
    sort!(matches)
    [m[1:end-length(after)] for m in matches]
end

"""
    config_to_toml(cfg::ExperimentConfig) -> String

Serialize an `ExperimentConfig` to a TOML-formatted string suitable for
storage in the `experiments.config` column or writing to disk. Uses the
stdlib `TOML.print` to handle quoting and escaping correctly.
"""
function config_to_toml(cfg::ExperimentConfig)::String
    function col_value(v)
        v isa Integer ? Int(v) : String(v)
    end
    # Omit nullable beamline fields when unset so a round-trip preserves
    # `nothing` instead of silently collapsing to 0.0.
    beamline = Dict{String,Any}()
    cfg.energy_kev    !== nothing && (beamline["energy_kev"]    = cfg.energy_kev)
    cfg.flight_path_m !== nothing && (beamline["flight_path_m"] = cfg.flight_path_m)
    d = Dict(
        "experiment" => Dict(
            "name"        => cfg.name,
            "description" => cfg.description,
            "manifest"    => cfg.manifest_file,
        ),
        "beamline" => beamline,
        "manifest" => Dict(
            "delimiter"      => cfg.delimiter,
            "skip_rows"      => cfg.skip_rows,
            "header_row"     => cfg.header_row,
            "sample_id"      => col_value(cfg.col_sample_id),
            "label"          => col_value(cfg.col_label),
            "name"           => col_value(cfg.col_name),
            "filenames"      => col_value(cfg.col_filenames),
            "notes_sample"   => col_value(cfg.col_notes_sample),
            "notes_exposure" => col_value(cfg.col_notes_exposure),
        ),
        "layout" => Dict(
            "data_dir"      => cfg.data_dir,
            "analysis_dir"  => cfg.analysis_dir,
            "exposure_type" => cfg.exposure_type,
        ),
        "files" => Dict(
            "integration" => cfg.integration_pattern,
            "image"       => cfg.image_pattern,
        ),
    )
    io = IOBuffer()
    TOML.print(io, d)
    String(take!(io))
end

"""
    config_from_db(db, experiment_id) -> ExperimentConfig

Read the stored TOML blob from `experiments.config` and parse it back into
an `ExperimentConfig`. If `config` is `NULL` (legacy experiment), falls back
to the built-in `simple` template, preserving backward compatibility.
"""
function config_from_db(db::SQLite.DB, experiment_id::Int)::ExperimentConfig
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT config FROM experiments WHERE id = ?", [experiment_id]))
    isempty(rows) && error("Experiment $experiment_id not found")
    blob = rows[1].config
    if blob === nothing || blob === missing
        return load_builtin_config("simple")
    end
    _build_config(TOML.parse(String(blob)))
end
