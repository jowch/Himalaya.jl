using HTTP, JSON3, DBInterface, Tables, Oxygen

_setup_infos(dir::AbstractString)::Vector{String} =
    isdir(dir) ?
        String[joinpath(dir, f) for f in (try readdir(dir) catch; String[] end)
               if startswith(f, "setup_info_") && endswith(f, ".txt")] :
        String[]

_listdir(dir::AbstractString) = try readdir(dir) catch; String[] end

_subdirs(d::AbstractString) = isdir(d) ?
    String[joinpath(d, x) for x in _listdir(d) if !startswith(x, ".") && isdir(joinpath(d, x))] :
    String[]

const LAYOUT_PROBE_STEMS = 9

function _layout_for_stem(stem::AbstractString, candidates::Vector{String})
    m = match(r"_0_\d+$", stem)
    suffix = m === nothing ? "" : m.match
    base   = m === nothing ? stem : stem[1:end - length(suffix)]
    for d in candidates
        if isfile(joinpath(d, base * "_tot.dat"))
            return (; dir = d,
                      image_pattern       = "{name}$(suffix).tif",
                      metadata_pattern    = "{name}$(suffix).prp",
                      integration_pattern = "{name}_tot.dat")
        end
    end
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

_is_tif(f) = endswith(f, ".tif") && !startswith(f, ".")

function _structural_ssrl_layout(stems, analysis_root)
    tot = joinpath(analysis_root, "automatic_analysis", "tot_files")
    isdir(tot) || return nothing
    suffixes = String[m.match for m in (match(r"_0_\d+$", s) for s in stems) if m !== nothing]
    isempty(suffixes) && return nothing
    counts = Dict{String,Int}()
    for x in suffixes; counts[x] = get(counts, x, 0) + 1; end
    suffix = first(sort!(collect(keys(counts)); by = x -> (-counts[x], x)))
    return (; dir = tot,
              image_pattern       = "{name}$(suffix).tif",
              metadata_pattern    = "{name}$(suffix).prp",
              integration_pattern = "{name}_tot.dat")
end

function _detect_integration_layout(data_dir::AbstractString, analysis_root::AbstractString)
    isdir(analysis_root) || return nothing
    tifs = sort!([f for f in _listdir(data_dir) if _is_tif(f)])
    isempty(tifs) && return nothing
    n = length(tifs)
    picks = unique(round.(Int, range(1, n; length = min(n, LAYOUT_PROBE_STEMS))))
    stems = [tifs[i][1:end-length(".tif")] for i in picks]
    candidates = String[analysis_root]
    append!(candidates, _subdirs(analysis_root))
    append!(candidates, _subdirs(joinpath(analysis_root, "automatic_analysis")))
    votes = Dict{Any,Int}()
    order = Any[]
    for stem in stems
        layout = _layout_for_stem(stem, candidates)
        layout === nothing && continue
        haskey(votes, layout) || push!(order, layout)
        votes[layout] = get(votes, layout, 0) + 1
    end
    isempty(order) && return _structural_ssrl_layout(stems, analysis_root)
    rank(l) = (votes[l], l.integration_pattern == "{name}_tot.dat")
    best = maximum(rank, order)
    return first(l for l in order if rank(l) == best)
end

"""
    resolve_experiment_layout(root) -> NamedTuple

Structural defaults for one setup folder: `(; name, data_dir, analysis_dir,
setup_file, setup_ambiguous, image_pattern, metadata_pattern, integration_pattern)`.
Only shallow reads, because the folder is usually on SMB.
"""
function resolve_experiment_layout(root::AbstractString)
    root = String(root)
    name = basename(rstrip(root, '/'))

    data_candidate = joinpath(root, "data")
    data_dir = isdir(data_candidate) ? data_candidate : root

    analysis_candidate = joinpath(root, "analysis")
    detected = _detect_integration_layout(data_dir, analysis_candidate)
    analysis_dir = detected !== nothing ? detected.dir :
                   (isdir(analysis_candidate) ? analysis_candidate : nothing)
    image_pattern       = detected !== nothing ? detected.image_pattern       : nothing
    metadata_pattern    = detected !== nothing ? detected.metadata_pattern    : nothing
    integration_pattern = detected !== nothing ? detected.integration_pattern : nothing

    setup_files = sort!(vcat(_setup_infos(root), _setup_infos(analysis_candidate)))
    setup_file = isempty(setup_files) ? nothing : last(setup_files)
    setup_ambiguous = isempty(setup_files) || length(unique(dirname.(setup_files))) > 1

    return (; name, data_dir, analysis_dir, setup_file, setup_ambiguous,
              image_pattern, metadata_pattern, integration_pattern)
end

function _pattern_stems(filelist, pat)
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

_is_near_ext(lext, patext) = lext == patext ||
    (startswith(lext, patext) && length(lext) == length(patext) + 1)

function _near_sidecar(filelist, pat, stem)
    occursin("{name}", pat) || return nothing
    pre, post = split(pat, "{name}"; limit = 2)
    pdot = findlast('.', post)
    pdot === nothing && return nothing
    patext = lowercase(post[pdot:end])
    target = pre * stem * post[1:prevind(post, pdot)]
    for f in filelist
        (startswith(f, target) && length(f) > length(target)) || continue
        _is_near_ext(lowercase(f[nextind(f, lastindex(target)):end]), patext) && return f
    end
    nothing
end

function _near_images(files, image_pattern)
    out = Dict{String,String}[]
    occursin("{name}", image_pattern) || return out
    ipre, ipost = split(image_pattern, "{name}"; limit = 2)
    pdot = findlast('.', ipost)
    pdot === nothing && return out
    patext = lowercase(ipost[pdot:end])
    patmid = ipost[1:prevind(ipost, pdot)]
    for f in files
        (startswith(f, ipre) && !endswith(f, ipost)) || continue
        rest = f[nextind(f, lastindex(ipre)):end]
        fdot = findlast('.', rest)
        fdot === nothing && continue
        fmid = rest[1:prevind(rest, fdot)]
        _is_near_ext(lowercase(rest[fdot:end]), patext) || continue
        endswith(fmid, patmid) || continue
        stem = chop(fmid; tail = length(patmid))
        isempty(stem) && continue
        push!(out, Dict("file" => stem, "miss" => "image", "near" => f))
    end
    out
end

const _GEOMETRY_MISSING = (energy_kev = nothing, energy_kev_source = "default",
    flight_path_m = nothing, flight_path_m_source = "default",
    beam_center_x = nothing, beam_center_x_source = "default",
    beam_center_y = nothing, beam_center_y_source = "default",
    pixel_size_um = nothing, pixel_size_um_source = "default")

const GEOMETRY_PRP_SAMPLE = 50

function _preview_geometry(data_dir, metadata_pattern, meta_stems, setup_file)
    prp_paths = if occursin("{name}", metadata_pattern)
        pre, post = split(metadata_pattern, "{name}"; limit = 2)
        sample = collect(Iterators.take(sort!(collect(meta_stems)), GEOMETRY_PRP_SAMPLE))
        String[joinpath(data_dir, pre * s * post) for s in sample]
    else
        String[]
    end
    setup_files = if setup_file !== nothing && isfile(setup_file)
        [setup_file]
    else
        ad = let a = joinpath(data_dir, "analysis"); isdir(a) ? a : data_dir end
        _setup_infos(ad)
    end
    geo, disc = try
        derive_geometry(prp_paths, setup_files)
    catch e
        @warn "preview: derive_geometry failed" data_dir exception = e
        (_GEOMETRY_MISSING, [])
    end
    _nm(x) = ismissing(x) ? nothing : x
    gd = Dict(
        :energy_kev            => _nm(geo.energy_kev),
        :energy_kev_source     => geo.energy_kev_source,
        :flight_path_m         => _nm(geo.flight_path_m),
        :flight_path_m_source  => geo.flight_path_m_source,
        :beam_center_x         => _nm(geo.beam_center_x),
        :beam_center_x_source  => geo.beam_center_x_source,
        :beam_center_y         => _nm(geo.beam_center_y),
        :beam_center_y_source  => geo.beam_center_y_source,
        :pixel_size_um         => _nm(geo.pixel_size_um),
        :pixel_size_um_source  => geo.pixel_size_um_source,
    )
    (gd, [Dict(:field => d.field, :message => d.message) for d in disc])
end

"""
    fs_manifest(data_dir; image_pattern, metadata_pattern, integration_pattern,
                analysis_dir, setup_file) -> NamedTuple

Read-only match of the three patterns against `data_dir` (integration against
`analysis_dir` when given) plus the geometry the scan would derive.
"""
function fs_manifest(data_dir::AbstractString;
        image_pattern::AbstractString = "{name}.tif",
        metadata_pattern::AbstractString = "{name}.prp",
        integration_pattern::AbstractString = "{name}.dat",
        analysis_dir::Union{AbstractString,Nothing} = nothing,
        setup_file::Union{AbstractString,Nothing} = nothing)
    files = filter(!startswith("."), _listdir(data_dir))
    integ_files = (analysis_dir !== nothing && isdir(analysis_dir)) ?
        filter(!startswith("."), _listdir(analysis_dir)) : files
    img   = _pattern_stems(files, image_pattern)
    meta  = _pattern_stems(files, metadata_pattern)
    integ = _pattern_stems(integ_files, integration_pattern)
    unmatched = Dict{String,String}[]
    for s in img,
        (label, set, flist, pat) in (("metadata", meta, files, metadata_pattern),
                                     ("integration", integ, integ_files, integration_pattern))
        s in set && continue
        near = _near_sidecar(flist, pat, s)
        push!(unmatched, near === nothing ?
            Dict("file" => s, "miss" => label) :
            Dict("file" => s, "miss" => label, "near" => near))
    end
    append!(unmatched, _near_images(files, image_pattern))
    geometry, discrepancies = _preview_geometry(data_dir, metadata_pattern, intersect(meta, img), setup_file)
    img_pre, img_post = occursin("{name}", image_pattern) ?
        split(image_pattern, "{name}"; limit = 2) : ("", "")
    names = sort!([img_pre * s * img_post for s in img])
    return (; total = length(files),
              image = length(img),
              metadata = length(meta),
              integration = length(integ),
              without_metadata = length(setdiff(img, meta)),
              without_trace = length(setdiff(img, integ)),
              unmatched, geometry, discrepancies,
              matched_files = names[1:min(12, length(names))])
end

_has_tifs(dir::AbstractString) = any(_is_tif, _listdir(dir))

"""
    setup_roots(path) -> Vector{String}

The setup folders `path` holds: `path` itself when its data folder has images,
else each immediate subfolder that does (the `<date>/setup1`, `setup2`, … layout).
"""
_setup_data_dir(root) = isdir(joinpath(root, "data")) ? joinpath(root, "data") : String(root)

function setup_roots(path::AbstractString)
    isdir(joinpath(path, "data")) && _has_tifs(joinpath(path, "data")) && return [String(path)]
    subs = String[d for d in sort!(_subdirs(path)) if _has_tifs(_setup_data_dir(d))]
    !isempty(subs) && return subs
    _has_tifs(path) ? [String(path)] : String[]
end

_canonical_dir(p) = (q = rstrip(String(p), '/'); isempty(q) ? "/" : (ispath(q) ? realpath(q) : q))

function experiment_using_dir(db, dir)
    plain = rstrip(String(dir), '/')
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, name FROM experiments WHERE rtrim(data_dir, '/') IN (?, ?) ORDER BY id LIMIT 1",
        [plain, _canonical_dir(dir)]))
    isempty(rows) ? nothing : (; id = Int(first(rows).id), name = String(first(rows).name))
end

_trim_slash(p) = rstrip(p, '/')

function _setup_name(path, root)
    parent = _trim_slash(path) == _trim_slash(root) ? dirname(_trim_slash(root)) : _trim_slash(path)
    pname = basename(parent)
    isempty(pname) ? basename(_trim_slash(root)) : "$(pname)/$(basename(_trim_slash(root)))"
end

function _setup_preview(db, path, root)
    r = resolve_experiment_layout(root)
    pats = (image       = something(r.image_pattern, "{name}.tif"),
            metadata    = something(r.metadata_pattern, "{name}.prp"),
            integration = something(r.integration_pattern, "{name}.dat"))
    m = fs_manifest(r.data_dir; image_pattern = pats.image, metadata_pattern = pats.metadata,
        integration_pattern = pats.integration, analysis_dir = r.analysis_dir,
        setup_file = r.setup_file)
    dup = experiment_using_dir(db, r.data_dir)
    Dict(:root => root,
         :name => _setup_name(path, root),
         :data_dir => r.data_dir,
         :analysis_dir => r.analysis_dir,
         :setup_file => r.setup_file,
         :setup_ambiguous => r.setup_ambiguous,
         :layout_detected => r.integration_pattern !== nothing,
         :patterns => Dict(:image => pats.image, :metadata => pats.metadata,
                           :integration => pats.integration),
         :exposures => m.image,
         :without_metadata => m.without_metadata,
         :without_trace => m.without_trace,
         :geometry => m.geometry,
         :discrepancies => m.discrepancies,
         :duplicate_of => dup === nothing ? nothing : Dict(:id => dup.id, :name => dup.name))
end

"""
    register_fs_routes!()

Read-only filesystem probes for the New-experiment screen: `suggest`, `manifest`,
and `preview` (every setup under a picked folder).
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
    @get "/api/fs/preview" function(req::HTTP.Request)
        q    = HTTP.queryparams(HTTP.URI(req.target))
        path = get(q, "path", "")
        isdir(path) || return _json(400, Dict(:error => "path is not a directory"))
        db = current_db()
        _json(200, Dict(:path => path,
                        :setups => [_setup_preview(db, path, root) for root in setup_roots(path)]))
    end
    @get "/api/fs/manifest" function(req::HTTP.Request)
        q    = HTTP.queryparams(HTTP.URI(req.target))
        path = get(q, "path", "")
        isdir(path) || return _json(400, Dict(:error => "path is not a directory"))
        opt(k) = (v = get(q, k, ""); isempty(v) ? nothing : v)
        m = fs_manifest(path;
            image_pattern       = get(q, "image_pattern", "{name}.tif"),
            metadata_pattern    = get(q, "metadata_pattern", "{name}.prp"),
            integration_pattern = get(q, "integration_pattern", "{name}.dat"),
            analysis_dir = opt("analysis_dir"),
            setup_file   = opt("setup_file"))
        _json(200, Dict(:total => m.total,
                   :matched => Dict(:image => m.image, :metadata => m.metadata, :integration => m.integration),
                   :unmatched => m.unmatched,
                   :geometry => m.geometry,
                   :discrepancies => m.discrepancies,
                   :matched_files => m.matched_files))
    end
end
