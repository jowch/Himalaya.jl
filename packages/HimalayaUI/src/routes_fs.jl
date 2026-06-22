using HTTP, JSON3, DBInterface, Tables, Oxygen

# ─────────────────────────────────────────────────────────────────────────────
# Filesystem probes for the ingest funnel (spec §6.1). Read-only — no DB
# writes, no with_idempotency, no apply_event!, no SSE, no client_op_id,
# no pendingDeferreds participation. Queue-orthogonal by design.
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# Structural experiment-layout resolver (funnel resolution, 2026-06-21).
# Given an experiment ROOT directory, auto-discover sensible defaults the user
# then corrects in Configuration before Approve. PURELY STRUCTURAL — no
# experiment.toml / manifest.toml (deprecated going forward); everything is
# inferred from the directory tree:
#   name         = basename(root)
#   data_dir     = the root or an immediate subdir holding the most .tif (raw images)
#   analysis_dir = the dir (under a non-data subtree) holding the most .dat (integration)
#   setup_file   = setup_info_*.txt found by a bounded tree walk (geometry source);
#                  it lives under analysis/, SEPARATE from the integration analysis_dir.
# Bounded for SMB latency: the data_dir subtree (the big one) is never walked;
# only sibling subdirs are walked, depth-capped.
# ─────────────────────────────────────────────────────────────────────────────

_setup_infos(dir::AbstractString)::Vector{String} =
    isdir(dir) ?
        String[joinpath(dir, f) for f in (try readdir(dir) catch; String[] end)
               if startswith(f, "setup_info_") && endswith(f, ".txt")] :
        String[]

# Bounded discovery of WHERE integration sidecars live AND which pattern names
# them. SSRL data nests integration under an analysis subtree, and the canonical
# integration is the per-base TOTAL: `analysis/automatic_analysis/tot_files/{base}_tot.dat`
# (`{base}` = the sample stem with its `_0_NNN` frame suffix stripped). So the
# scan's flat `sidecar(analysis_dir, stem, dat_pattern)` only finds them when
# analysis_dir points at the real `.dat` dir AND the patterns match the convention.
#
# Probe one sample data stem against shallow candidate dirs (analysis, its
# immediate subdirs, and the immediate subdirs of `analysis/automatic_analysis`):
#   1. TOT (canonical): `{base}_tot.dat` → patterns key off the frame suffix.
#   2. Per-frame `.dat` (e.g. sastool): `{stem}.dat` → default `{name}` patterns.
# Returns `(; dir, image_pattern, metadata_pattern, integration_pattern)` or
# `nothing` (caller falls back to root/analysis with default patterns).
# SMB-bounded: one `readdir(data_dir)` (as validate/manifest already do) for a
# sample stem, then a handful of `isdir`/`isfile` stats — NO deep walk.
function _detect_integration_layout(data_dir::AbstractString, analysis_root::AbstractString)
    isdir(analysis_root) || return nothing
    sample_tif = nothing
    for f in (try readdir(data_dir) catch; String[] end)
        if endswith(f, ".tif"); sample_tif = f; break; end
    end
    sample_tif === nothing && return nothing
    stem = sample_tif[1:end-length(".tif")]
    # Frame suffix (e.g. "_0_001") if present; `{base}` is the stem without it.
    m = match(r"_0_\d+$", stem)
    suffix = m === nothing ? "" : m.match
    base   = m === nothing ? stem : stem[1:end - length(suffix)]

    subdirs(d) = isdir(d) ?
        String[joinpath(d, x) for x in (try readdir(d) catch; String[] end)
               if isdir(joinpath(d, x))] : String[]
    candidates = String[analysis_root]
    append!(candidates, subdirs(analysis_root))
    append!(candidates, subdirs(joinpath(analysis_root, "automatic_analysis")))

    # 1. Canonical TOTAL integration: {base}_tot.dat. The image/metadata patterns
    #    key off the actual frame suffix so the scan's stem resolves to {base}.
    for d in candidates
        if isfile(joinpath(d, base * "_tot.dat"))
            return (; dir = d,
                      image_pattern       = "{name}$(suffix).tif",
                      metadata_pattern    = "{name}$(suffix).prp",
                      integration_pattern = "{name}_tot.dat")
        end
    end
    # 2. Per-frame integration: {stem}.dat (default {name} patterns).
    for d in candidates
        if isfile(joinpath(d, stem * ".dat"))
            return (; dir = d,
                      image_pattern       = "{name}.tif",
                      metadata_pattern    = "{name}.prp",
                      integration_pattern = "{name}.dat")
        end
    end
    return nothing
end

"""
    resolve_experiment_layout(root) -> NamedTuple

Structural auto-discovery for the ingest funnel. Returns
`(; name, data_dir, analysis_dir, setup_file, setup_ambiguous)`. All values are
DEFAULTS the user corrects in Configuration — so this stays CHEAP (a handful of
shallow `isdir`/`readdir` calls, NEVER a deep walk and NEVER a full read of the
huge raw-image directory), because it runs interactively over slow SMB.

Heuristics (name-based, with a correctable safety net):
- `name`         = `basename(root)`.
- `data_dir`     = `root/data` if it exists, else `root` (raw images usually
                   live in a `data` subdir; the user corrects if not).
- `analysis_dir` = the directory that actually holds the integration `.dat`
                   sidecars, discovered by probing a sample stem against the
                   analysis subtree (e.g. `analysis/automatic_analysis/sastool`);
                   falls back to `root/analysis`, else `nothing`. The user can
                   still correct it in Configuration.
- `setup_file`   = the latest `setup_info_*.txt` found at the TOP of `root` or
                   `root/analysis` (where it actually lives — NOT next to the
                   integration files). This is the geometry source, so it is
                   resolved precisely; `setup_ambiguous` flags none/multiple.
"""
function resolve_experiment_layout(root::AbstractString)
    root = String(root)
    name = basename(rstrip(root, '/'))

    data_candidate = joinpath(root, "data")
    data_dir = isdir(data_candidate) ? data_candidate : root

    analysis_candidate = joinpath(root, "analysis")
    # Discover the directory that actually holds the integration `.dat` files
    # (often nested, e.g. analysis/automatic_analysis/tot_files) AND the file
    # patterns that name them, so the scan finds them. Fall back to root/analysis
    # with default `{name}` patterns when nothing is discovered.
    detected = _detect_integration_layout(data_dir, analysis_candidate)
    analysis_dir = detected !== nothing ? detected.dir :
                   (isdir(analysis_candidate) ? analysis_candidate : nothing)
    image_pattern       = detected !== nothing ? detected.image_pattern       : nothing
    metadata_pattern    = detected !== nothing ? detected.metadata_pattern    : nothing
    integration_pattern = detected !== nothing ? detected.integration_pattern : nothing

    # setup_info lives at the top of root or root/analysis — two shallow reads.
    setup_files = sort!(vcat(_setup_infos(root), _setup_infos(analysis_candidate)))
    setup_file = isempty(setup_files) ? nothing : last(setup_files)  # latest by name
    setup_ambiguous = isempty(setup_files) || length(unique(dirname.(setup_files))) > 1

    return (; name, data_dir, analysis_dir, setup_file, setup_ambiguous,
              image_pattern, metadata_pattern, integration_pattern)
end

"""
    register_fs_routes!()

Filesystem probes for the ingest funnel (spec §6.1). Read-only, no DB writes.
- GET /api/fs/suggest?prefix=  → directory autocomplete for the picker.
- GET /api/fs/validate?path=   → cheap picker gate (exists + not-already-an-experiment).
- GET /api/fs/manifest?path=&{image,metadata,integration}_pattern=  → phase-1 file manifest.

`validate` and `manifest` are added in Tasks 1.4 / 1.5.
"""
function register_fs_routes!()
    @get "/api/fs/suggest" function(req::HTTP.Request)
        q      = HTTP.queryparams(HTTP.URI(req.target))
        prefix = get(q, "prefix", "")
        isempty(prefix) && return _json(200, Dict(:suggestions => String[]))
        dir  = isdir(prefix) ? prefix : dirname(prefix)
        base = isdir(prefix) ? "" : basename(prefix)
        isdir(dir) || return _json(200, Dict(:suggestions => String[]))
        kids = String[]
        for name in readdir(dir; sort = true)
            startswith(name, base) || continue
            full = joinpath(dir, name)
            isdir(full) && push!(kids, full)
            length(kids) >= 20 && break
        end
        _json(200, Dict(:suggestions => kids))
    end
    @get "/api/fs/resolve" function(req::HTTP.Request)
        q    = HTTP.queryparams(HTTP.URI(req.target))
        path = get(q, "path", "")
        isdir(path) || return _json(400, Dict(:error => "path is not a directory"))
        r = resolve_experiment_layout(path)
        _json(200, Dict(:name => r.name,
                   :data_dir => r.data_dir,
                   :analysis_dir => r.analysis_dir,
                   :setup_file => r.setup_file,
                   :setup_ambiguous => r.setup_ambiguous,
                   :image_pattern => r.image_pattern,
                   :metadata_pattern => r.metadata_pattern,
                   :integration_pattern => r.integration_pattern))
    end
    @get "/api/fs/validate" function(req::HTTP.Request)
        q    = HTTP.queryparams(HTTP.URI(req.target))
        path = get(q, "path", "")
        if isempty(path) || !isdir(path)
            return _json(200, Dict(:ok => false, :matched => 0, :scanned => 0,
                               :message => "path does not exist or is not a directory"))
        end
        dup = !isempty(DBInterface.execute(current_db(),
            "SELECT 1 FROM experiments WHERE data_dir = ? LIMIT 1", [path]) |> Tables.rowtable)
        if dup
            return _json(200, Dict(:ok => false, :matched => 0, :scanned => 0,
                               :message => "an experiment already uses this directory"))
        end
        scanned = count(!startswith("."), readdir(path))   # cheap; rich count is /manifest
        _json(200, Dict(:ok => true, :matched => scanned, :scanned => scanned, :message => nothing))
    end
    @get "/api/fs/manifest" function(req::HTTP.Request)
        q    = HTTP.queryparams(HTTP.URI(req.target))
        path = get(q, "path", "")
        isdir(path) || return _json(400, Dict(:error => "path is not a directory"))
        pats = (image       = get(q, "image_pattern", "{name}.tif"),
                metadata    = get(q, "metadata_pattern", "{name}.prp"),
                integration = get(q, "integration_pattern", "{name}.dat"))
        files = filter(!startswith("."), readdir(path))
        # Integration (.dat) lives in the ANALYSIS subtree, not data_dir — the real
        # scan pairs it via `sidecar(analysis_dir, stem, dat_pattern)` (grouping.jl
        # `scan_directory`). Mirror that here so the preview's integration count is
        # truthful: match `integration_pattern` against `analysis_dir` when supplied.
        # Absent → fall back to data_dir (backward-compatible with older callers/tests).
        analysis_dir_q = get(q, "analysis_dir", "")
        integ_files = (!isempty(analysis_dir_q) && isdir(analysis_dir_q)) ?
            filter(!startswith("."), readdir(analysis_dir_q)) : files
        # {name}-capture per type, inlined. Do NOT add a module helper — grouping.jl
        # / config.jl already own pattern matching; if the manifest's needs grow,
        # call the shared matcher (config.jl `_matches_prefix_with_boundary` /
        # `resolve_files`) rather than maintaining a second one (ponytail review).
        function stems(filelist, pat)
            occursin("{name}", pat) || return Set{String}()
            pre, post = split(pat, "{name}"; limit = 2)
            out = Set{String}()
            for f in filelist
                (startswith(f, pre) && endswith(f, post) &&
                    length(f) > length(pre) + length(post)) || continue
                push!(out, f[nextind(f, lastindex(pre)):prevind(f, lastindex(f) - length(post) + 1)])
            end
            out
        end
        img  = stems(files, pats.image)
        meta = stems(files, pats.metadata)
        integ = stems(integ_files, pats.integration)
        unmatched = Dict{String,String}[]
        for s in img, (label, set) in (("metadata", meta), ("integration", integ))
            s in set || push!(unmatched, Dict("file" => s, "miss" => label))
        end

        # ── Geometry preview (read-only) ──
        # `path` is the data_dir (the frontend passes the resolved/edited data_dir).
        # `setup_file` (optional) is the geometry source resolved by /api/fs/resolve;
        # when supplied we use it directly (it lives under analysis/, NOT next to the
        # integration files). Absent → fall back to the old data_dir/analysis scan so
        # older callers/tests keep working.
        data_dir       = path
        setup_file_q   = get(q, "setup_file", "")

        # Reconstruct matched metadata filenames from stems (same pre/post split).
        _nm(x) = ismissing(x) ? nothing : x
        geo_missing = (energy_kev = nothing, energy_kev_source = "default",
                       flight_path_m = nothing, flight_path_m_source = "default",
                       beam_center_x = nothing, beam_center_x_source = "default",
                       beam_center_y = nothing, beam_center_y_source = "default",
                       pixel_size_um = nothing, pixel_size_um_source = "default")

        geo_dict, disc_list = let
            prp_paths = if occursin("{name}", pats.metadata)
                meta_pre, meta_post = split(pats.metadata, "{name}"; limit = 2)
                # Geometry is CONSTANT per run, so sample up to 50 PRPs — keeps the
                # preview fast over SMB (a full run is thousands of files, and each
                # would be an isfile + read round-trip). The first 50 suffice.
                sample = collect(Iterators.take(meta, 50))
                [joinpath(data_dir, meta_pre * s * meta_post) for s in sample
                    if isfile(joinpath(data_dir, meta_pre * s * meta_post))]
            else
                String[]
            end
            setup_files = if !isempty(setup_file_q) && isfile(setup_file_q)
                [setup_file_q]
            else
                # Fallback (no explicit setup_file): old data_dir/analysis scan.
                ad = let a = joinpath(data_dir, "analysis"); isdir(a) ? a : data_dir end
                isdir(ad) ?
                    [joinpath(ad, f) for f in readdir(ad)
                        if startswith(f, "setup_info_") && endswith(f, ".txt")] :
                    String[]
            end
            # ponytail: if prp_paths grows to thousands, consider sampling ~100
            # representative files for geometry preview (geometry is constant per run).
            geo_result, disc_result = try
                derive_geometry(prp_paths, setup_files)
            catch e
                @warn "manifest: derive_geometry failed" path exception=e
                (geo_missing, [])
            end
            gd = Dict(
                :energy_kev            => _nm(geo_result.energy_kev),
                :energy_kev_source     => geo_result.energy_kev_source,
                :flight_path_m         => _nm(geo_result.flight_path_m),
                :flight_path_m_source  => geo_result.flight_path_m_source,
                :beam_center_x         => _nm(geo_result.beam_center_x),
                :beam_center_x_source  => geo_result.beam_center_x_source,
                :beam_center_y         => _nm(geo_result.beam_center_y),
                :beam_center_y_source  => geo_result.beam_center_y_source,
                :pixel_size_um         => _nm(geo_result.pixel_size_um),
                :pixel_size_um_source  => geo_result.pixel_size_um_source,
            )
            dd = [Dict(:field => d.field, :message => d.message) for d in disc_result]
            (gd, dd)
        end

        # ── matched_files: capped sample of image filenames (for "latest files" card) ──
        matched_files = let
            img_pre, img_post = occursin("{name}", pats.image) ?
                split(pats.image, "{name}"; limit = 2) : ("", "")
            names = sort!([img_pre * s * img_post for s in img])
            names[1:min(12, length(names))]
        end

        _json(200, Dict(:total => length(files),
                   :matched => Dict(:image => length(img), :metadata => length(meta), :integration => length(integ)),
                   :unmatched => unmatched,
                   :geometry => geo_dict,
                   :discrepancies => disc_list,
                   :matched_files => matched_files))
    end
end
