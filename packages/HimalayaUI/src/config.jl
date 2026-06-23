using TOML

const VALID_EXPOSURE_TYPES = ("simple",)

"""
    ExperimentConfig

Parsed representation of an `experiment.toml` file. Drives all file I/O conventions
for an experiment: manifest column layout, directory paths, file name patterns,
beamline parameters.

Supported TOML sections:
- `[experiment]`: name, description, manifest (relative path to manifest CSV)
- `[beamline]`: energy_kev, flight_path_m, q_units, beam_center_x, beam_center_y, pixel_size_um
- `[manifest]`: delimiter, skip_rows, header_row, sample_id, name, display_name,
  filenames, notes_sample, notes_exposure (each column = Int index or String header name).
  `header_row = 0` is the sentinel for "no header row; columns are positional".
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
    q_units            ::String
    beam_center_x      ::Union{Float64,Nothing}
    beam_center_y      ::Union{Float64,Nothing}
    pixel_size_um      ::Union{Float64,Nothing}
    # [manifest]
    delimiter          ::String
    skip_rows          ::Int
    header_row         ::Int
    col_sample_id      ::Union{Int,String}
    col_name           ::Union{Int,String}
    col_display_name   ::Union{Int,String}
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

    if haskey(mf, "label")
        error("deprecated key '[manifest].label' in experiment.toml. " *
              "The manifest column meanings were swapped: column 2 is now `name` " *
              "(stable identifier), column 3 is now `display_name` (user-friendly label). " *
              "Run `himalaya migrate-toml <experiment-dir>` to upgrade automatically.")
    end

    ExperimentConfig(
        get(exp, "name",        ""),
        get(exp, "description", ""),
        get(exp, "manifest",    "manifest.csv"),
        get(bl,  "energy_kev",    nothing),
        get(bl,  "flight_path_m", nothing),
        get(bl,  "q_units",       "A-1"),
        get(bl,  "beam_center_x", nothing),
        get(bl,  "beam_center_y", nothing),
        get(bl,  "pixel_size_um", nothing),
        get(mf,  "delimiter",      "\t"),
        get(mf,  "skip_rows",      1),
        get(mf,  "header_row",     0),
        _coerce_col(get(mf, "sample_id",      1),  "manifest.sample_id"),
        _coerce_col(get(mf, "name",         2),  "manifest.name"),
        _coerce_col(get(mf, "display_name", 3),  "manifest.display_name"),
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
    _matches_prefix_with_boundary(f, prefix, suffix) -> Bool

Like `startswith(f, prefix) && endswith(f, suffix)`, but rejects matches where
the character immediately after `prefix` is alphanumeric. This prevents
prefix collisions where one manifest entry is a string-prefix of another:
`JC_C04` would otherwise catch every `JC_C04s_*.dat` file too. With the
boundary rule, the next character must be non-alphanumeric (typically `_`,
`.`, or `-`) — or the string must end exactly at the prefix.
"""
function _matches_prefix_with_boundary(f::AbstractString, prefix::AbstractString,
                                       suffix::AbstractString)::Bool
    startswith(f, prefix) || return false
    endswith(f, suffix)   || return false
    plen = ncodeunits(prefix)
    ncodeunits(f) == plen && return true   # exact match (suffix may be "")
    # Filenames are ASCII in practice (SAXS detector output), but guard against
    # mid-codepoint indexing if a non-ASCII filename ever slips in: `isvalid`
    # confirms `plen+1` lands on a char boundary before we read it.
    isvalid(f, plen + 1) || return false
    c = f[plen + 1]
    !(isletter(c) || isdigit(c))
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
        _matches_prefix_with_boundary(f, file_prefix_local, after)
    end
    sort!(matches)
    [m[1:end-length(after)] for m in matches]
end

"""
    resolve_file_path(cfg, base_dir, prefix, file_pattern) -> Union{String,Nothing}

Like [`resolve_files`](@ref) but returns the absolute path of the first match
(sorted), or `nothing` if no file matches. Used for the image-pattern lookup
where instrumentation typically writes one image per integration stem with
trailing tokens that aren't in the manifest (e.g. integration stem
`JC_D01_1_S2449` ↔ image `JC_D01_1_S2449_0_001.tiff`).
"""
function resolve_file_path(
    ::ExperimentConfig,
    base_dir::AbstractString,
    prefix::AbstractString,
    file_pattern::String,
)::Union{String,Nothing}
    parts = split(file_pattern, "{name}"; limit=2)
    length(parts) == 2 || error("file pattern must contain exactly one {name}: $file_pattern")
    before, after = String(parts[1]), String(parts[2])

    scan_subdir       = dirname(before)
    file_prefix_local = basename(before) * prefix
    scan_dir = isempty(scan_subdir) ? String(base_dir) : joinpath(base_dir, scan_subdir)

    isdir(scan_dir) || return nothing

    matches = filter(readdir(scan_dir)) do f
        _matches_prefix_with_boundary(f, file_prefix_local, after)
    end
    isempty(matches) && return nothing
    sort!(matches)
    joinpath(scan_dir, matches[1])
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
    cfg.beam_center_x !== nothing && (beamline["beam_center_x"] = cfg.beam_center_x)
    cfg.beam_center_y !== nothing && (beamline["beam_center_y"] = cfg.beam_center_y)
    cfg.pixel_size_um !== nothing && (beamline["pixel_size_um"] = cfg.pixel_size_um)
    beamline["q_units"] = cfg.q_units
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
            "name"           => col_value(cfg.col_name),
            "display_name"   => col_value(cfg.col_display_name),
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
    migrate_manifest_toml_text(text) -> (new_text, changed)

Pure-text rewrite of a TOML blob from the legacy `[manifest].label`/`name`
shape (column 2 is the stable identifier, column 3 is the user-facing label)
to the canonical `[manifest].name`/`display_name` shape (column 2 is the
stable identifier `name`, column 3 is the editable `display_name`). See
PR #107.

Section-aware: only substitutes inside the `[manifest]` block so a
`name` key in `[experiment]` (the human-readable experiment name) is not
touched. Per-line state machine: a single line is EITHER a `label=` rewrite
OR a `name=` rewrite, never both — so `name` doesn't catch the line we just
rewrote from `label`.

Idempotent: returns `(text, false)` when no rewrite is needed (already
migrated, or no `[manifest].label` present). Errors when the blob has
both `label` and `display_name` in `[manifest]` — that's a corrupt state
that needs operator attention; better to surface the corruption than to
guess.
"""
function migrate_manifest_toml_text(text::AbstractString)::Tuple{String,Bool}
    lines = collect(eachline(IOBuffer(text); keep=true))

    # First pass: classify [manifest] contents.
    section = ""
    has_label = false
    has_display_name = false
    for line in lines
        m = match(r"^\s*\[([A-Za-z0-9_]+)\]\s*$", line)
        if m !== nothing
            section = m.captures[1]; continue
        end
        if section == "manifest"
            occursin(r"^\s*label\s*=", line)        && (has_label = true)
            occursin(r"^\s*display_name\s*=", line) && (has_display_name = true)
        end
    end

    # No-op cases.
    has_label || return (String(text), false)
    if has_label && has_display_name
        error("manifest TOML has both `label` and `display_name` in [manifest]; manual edit needed")
    end

    # Second pass: rewrite only inside [manifest].
    section = ""
    out = IOBuffer()
    for line in lines
        m = match(r"^\s*\[([A-Za-z0-9_]+)\]\s*$", line)
        if m !== nothing
            section = m.captures[1]
            print(out, line); continue
        end
        if section == "manifest"
            if (m2 = match(r"^(\s*)label(\s*=\s*\S+)(.*)$", line)) !== nothing
                rest = rstrip(m2.captures[3], '\n')
                # Preserve trailing newline iff present in the original.
                nl = endswith(line, '\n') ? "\n" : ""
                print(out, m2.captures[1] * "name" * m2.captures[2] * rest * nl)
            elseif (m3 = match(r"^(\s*)name(\s*=\s*\S+)(.*)$", line)) !== nothing
                rest = rstrip(m3.captures[3], '\n')
                nl = endswith(line, '\n') ? "\n" : ""
                print(out, m3.captures[1] * "display_name" * m3.captures[2] * rest * nl)
            else
                print(out, line)
            end
        else
            print(out, line)
        end
    end
    (String(take!(out)), true)
end

"""
    config_from_db(db, experiment_id) -> ExperimentConfig

Read the stored TOML blob from `experiments.config` and parse it back into
an `ExperimentConfig`. If `config` is `NULL` (legacy experiment), falls back
to the built-in `simple` template, preserving backward compatibility.
"""
function config_from_db(db::SQLite.DB, experiment_id::Int)::ExperimentConfig
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT config, image_pattern, integration_pattern FROM experiments WHERE id = ?",
        [experiment_id]))
    isempty(rows) && error("Experiment $experiment_id not found")
    r = rows[1]
    base = (r.config === nothing || r.config === missing) ?
        load_builtin_config("simple") : _build_config(TOML.parse(String(r.config)))
    # The per-experiment pattern COLUMNS are the source of truth for HTTP-ingested
    # experiments: that path stores patterns in the columns and leaves the TOML
    # `config` blob NULL (the blob is deprecated). Without this override,
    # config_from_db returns the builtin `{name}.dat` and analyze_exposure! can
    # never resolve a real `_tot.dat` integration trace, so nothing indexes.
    _apply_db_patterns(base, r.image_pattern, r.integration_pattern)
end

"""
    _apply_db_patterns(cfg, image_col, integration_col) -> ExperimentConfig

Override a base config's `image_pattern` / `integration_pattern` with the
experiment's column values when they are present and well-formed (contain
`{name}`). The columns win over the deprecated TOML blob. Returns `cfg` unchanged
when both columns are absent. Robust to field reordering (rebuilds positionally
from `fieldnames`).
"""
function _apply_db_patterns(cfg::ExperimentConfig, image_col, integration_col)::ExperimentConfig
    pick(col, cur) = (col !== nothing && col !== missing &&
                      occursin("{name}", String(col))) ? String(col) : cur
    img   = pick(image_col, cfg.image_pattern)
    integ = pick(integration_col, cfg.integration_pattern)
    (img == cfg.image_pattern && integ == cfg.integration_pattern) && return cfg
    fn   = fieldnames(ExperimentConfig)
    vals = Any[getfield(cfg, f) for f in fn]
    vals[findfirst(==(:image_pattern), fn)]       = img
    vals[findfirst(==(:integration_pattern), fn)] = integ
    ExperimentConfig(vals...)
end
