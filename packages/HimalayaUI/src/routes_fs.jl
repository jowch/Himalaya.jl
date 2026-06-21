using HTTP, JSON3, DBInterface, Tables, Oxygen

# ─────────────────────────────────────────────────────────────────────────────
# Filesystem probes for the ingest funnel (spec §6.1). Read-only — no DB
# writes, no with_idempotency, no apply_event!, no SSE, no client_op_id,
# no pendingDeferreds participation. Queue-orthogonal by design.
# ─────────────────────────────────────────────────────────────────────────────

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
        # {name}-capture per type, inlined. Do NOT add a module helper — grouping.jl
        # / config.jl already own pattern matching; if the manifest's needs grow,
        # call the shared matcher (config.jl `_matches_prefix_with_boundary` /
        # `resolve_files`) rather than maintaining a second one (ponytail review).
        function stems(pat)
            occursin("{name}", pat) || return Set{String}()
            pre, post = split(pat, "{name}"; limit = 2)
            out = Set{String}()
            for f in files
                (startswith(f, pre) && endswith(f, post) &&
                    length(f) > length(pre) + length(post)) || continue
                push!(out, f[nextind(f, lastindex(pre)):prevind(f, lastindex(f) - length(post) + 1)])
            end
            out
        end
        img, meta, integ = stems(pats.image), stems(pats.metadata), stems(pats.integration)
        unmatched = Dict{String,String}[]
        for s in img, (label, set) in (("metadata", meta), ("integration", integ))
            s in set || push!(unmatched, Dict("file" => s, "miss" => label))
        end

        # ── Geometry preview (read-only; mirrors ingest.jl + routes_experiments.jl) ──
        data_dir     = path
        analysis_dir = let ad = joinpath(data_dir, "analysis"); isdir(ad) ? ad : data_dir end

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
                [joinpath(data_dir, meta_pre * s * meta_post) for s in meta
                    if isfile(joinpath(data_dir, meta_pre * s * meta_post))]
            else
                String[]
            end
            setup_files = isdir(analysis_dir) ?
                [joinpath(analysis_dir, f) for f in readdir(analysis_dir)
                    if startswith(f, "setup_info_") && endswith(f, ".txt")] :
                String[]
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
