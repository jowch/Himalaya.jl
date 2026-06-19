# packages/HimalayaUI/src/grouping.jl
#
# Directory scanning + sample grouping (spec §5, §9.1).
#
# scan_directory: enumerate (tif, prp, dat) triplets from a beamtime directory.
# group_into_samples: cluster exposures into loads and samples (Tasks 5–7).
# scan_and_group!: transactional orchestrator (ingest.jl, Task 8).
#
# Naming note: `auto_group` already exists in pipeline.jl for index-candidate
# grouping (an unrelated concept). This module uses `group_into_samples` (spec §9.1).

"""
    ExposureMeta

Per-file metadata collected during directory scan. All fields except `stem`
are nullable (`nothing` or `missing`).
"""
struct ExposureMeta
    stem         ::String
    tif_path     ::Union{String, Nothing}
    prp_path     ::Union{String, Nothing}
    dat_path     ::Union{String, Nothing}
    prp          ::Union{NamedTuple, Nothing}  # parse_prp result; nothing if no .prp
end

Base.haskey(m::ExposureMeta, k::Symbol) = k in fieldnames(ExposureMeta)

"""
    scan_directory(data_dir, analysis_dir;
                   tif_pattern = "{name}.tif",
                   prp_pattern = "{name}.prp",
                   dat_pattern = "{name}.dat") -> Vector{ExposureMeta}

Enumerate every TIFF found in `data_dir`, then pair each with its PRP and .dat
sidecar. Returns one `ExposureMeta` per TIFF stem, sorted by stem.

Reuses the `resolve_file_path` logic from `config.jl` for sidecar lookup once
a stem is known, constructing a minimal `ExperimentConfig` for dispatch (resolves
the §9.1 open question: loosen dispatch vs. construct minimal config — we construct
a minimal config since `resolve_file_path` takes `ExperimentConfig` for dispatch
only and none of its fields are read by the dispatch branch).

TIF enumeration does NOT use `resolve_files` with an empty prefix because
`_matches_prefix_with_boundary` rejects filenames whose first character is
alphanumeric when the prefix is empty — which is the common case. Instead we
scan the directory directly and strip the TIF suffix to derive stems.
"""
function scan_directory(
    data_dir::AbstractString,
    analysis_dir::AbstractString;
    tif_pattern::String = "{name}.tif",
    prp_pattern::String = "{name}.prp",
    dat_pattern::String = "{name}.dat",
)::Vector{ExposureMeta}
    # Build a minimal ExperimentConfig for the dispatch arg to resolve_file_path.
    # The dispatch method only needs the struct type; no fields are read.
    cfg = _minimal_scan_config(tif_pattern, prp_pattern, dat_pattern,
                               String(data_dir), String(analysis_dir))

    # Derive the literal TIF suffix from the pattern (e.g. "{name}.tif" → ".tif").
    tif_parts = split(tif_pattern, "{name}"; limit=2)
    length(tif_parts) == 2 || error("tif_pattern must contain {{name}}: $tif_pattern")
    tif_prefix_literal = String(tif_parts[1])  # e.g. "" for "{name}.tif"
    tif_suffix = String(tif_parts[2])           # e.g. ".tif"

    # Enumerate all files in data_dir that match the TIF pattern (exact suffix).
    scan_subdir = dirname(tif_prefix_literal)
    tif_scan_dir = isempty(scan_subdir) ? String(data_dir) : joinpath(data_dir, scan_subdir)

    if !isdir(tif_scan_dir)
        return ExposureMeta[]
    end

    tif_files = filter(readdir(tif_scan_dir)) do f
        startswith(f, tif_prefix_literal) && endswith(f, tif_suffix)
    end
    sort!(tif_files)

    # Strip the TIF suffix to get stems, then pair each stem with its sidecars.
    suffix_len = length(tif_suffix)
    metas = ExposureMeta[]
    for fname in tif_files
        stem = fname[1:end - suffix_len]
        tif_path  = joinpath(tif_scan_dir, fname)
        prp_path  = resolve_file_path(cfg, data_dir, stem, prp_pattern)
        dat_path  = resolve_file_path(cfg, analysis_dir, stem, dat_pattern)
        prp_parsed = prp_path !== nothing ? parse_prp(prp_path) : nothing
        push!(metas, ExposureMeta(stem, tif_path, prp_path, dat_path, prp_parsed))
    end
    return metas
end

"""
    _minimal_scan_config(tif_pattern, prp_pattern, dat_pattern, data_dir, analysis_dir) -> ExperimentConfig

Construct the minimal `ExperimentConfig` needed to dispatch `resolve_files` /
`resolve_file_path`. The manifest and beamline fields are set to safe sentinel
values and are never read during a scan.

Field order confirmed from config.jl struct definition (23 fields total):
  name, description, manifest_file                         # [experiment]
  energy_kev, flight_path_m, q_units,
    beam_center_x, beam_center_y, pixel_size_um           # [beamline]
  delimiter, skip_rows, header_row,
    col_sample_id, col_name, col_display_name,
    col_filenames, col_notes_sample, col_notes_exposure    # [manifest]
  data_dir, analysis_dir, exposure_type                    # [layout]
  integration_pattern, image_pattern                       # [files]

If ExperimentConfig gains new fields, this positional call must be updated
in lockstep — re-read config.jl before editing.
"""
function _minimal_scan_config(
    tif_pattern::String,
    prp_pattern::String,
    dat_pattern::String,
    data_dir::String,
    analysis_dir::String,
)::ExperimentConfig
    # Positional constructor — ExperimentConfig has no @kwdef / keyword form.
    ExperimentConfig(
        # [experiment]
        "", "", "",
        # [beamline]
        nothing, nothing, "A-1", nothing, nothing, nothing,
        # [manifest]
        ",", 0, 0, 1, 1, 1, 1, 1, 1,
        # [layout]
        data_dir, analysis_dir, "simple",
        # [files]  integration_pattern=dat, image_pattern=tif
        dat_pattern, tif_pattern,
    )
end
