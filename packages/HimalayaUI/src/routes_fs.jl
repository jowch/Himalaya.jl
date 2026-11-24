using HTTP, DBInterface, Tables, Oxygen

_setup_infos(dir::AbstractString)::Vector{String} =
    String[joinpath(dir, f) for f in list_names(dir) if startswith(f, "setup_info_") && endswith(f, ".txt")]

_listdir(dir::AbstractString) = list_names(dir)

_subdirs(d::AbstractString) = String[joinpath(d, x) for x in folder_names(d)]

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

function resolve_experiment_layout(root::AbstractString)
    root = String(root)
    name = basename(rstrip(root, '/'))

    data_candidate = joinpath(root, "data")
    data_dir = isdir(data_candidate) ? data_candidate : root

    analysis_candidate = joinpath(root, "analysis")
    mode = detect_setup_mode(data_dir)
    detected = mode == "autosampler" ? autosampler_layout(root) :
               _detect_integration_layout(data_dir, analysis_candidate)
    analysis_dir = detected !== nothing && detected.dir !== nothing ? detected.dir :
                   (isdir(analysis_candidate) ? analysis_candidate : nothing)
    image_pattern       = detected !== nothing ? detected.image_pattern       : nothing
    metadata_pattern    = detected !== nothing ? detected.metadata_pattern    : nothing
    integration_pattern = detected !== nothing ? detected.integration_pattern : nothing

    setup_files = sort!(vcat(_setup_infos(root), _setup_infos(analysis_candidate)))
    setup_file = isempty(setup_files) ? nothing : last(setup_files)
    setup_ambiguous = isempty(setup_files) || length(unique(dirname.(setup_files))) > 1

    return (; name, mode, data_dir, analysis_dir, setup_file, setup_ambiguous,
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

_has_tifs(dir::AbstractString) = any_name(_is_tif, dir)

_setup_data_dir(root) = isdir(joinpath(root, "data")) ? joinpath(root, "data") : String(root)

function setup_roots(path::AbstractString)
    isdir(joinpath(path, "data")) && _has_tifs(joinpath(path, "data")) && return [String(path)]
    subs = String[d for d in natural_sort(_subdirs(path)) if _has_tifs(_setup_data_dir(d))]
    !isempty(subs) && return subs
    _has_tifs(path) ? [String(path)] : String[]
end

_canonical_dir(p) = (q = rstrip(String(p), '/'); isempty(q) ? "/" : (ispath(q) ? realpath(q) : q))

_experiment_label(name, path) =
    !ismissing(name) && !isempty(name) ? String(name) : ismissing(path) ? "" : basename(_trim_slash(String(path)))

_named_experiment(row) = (; id = Int(row.id), name = _experiment_label(row.name, row.path))

function experiment_using_dir(db, dir)
    plain = rstrip(String(dir), '/')
    rows = sqlrows(DBInterface.execute(db,
        """SELECT x.id, x.name, x.path FROM setups st JOIN experiments x ON x.id = st.experiment_id
            WHERE rtrim(st.data_dir, '/') IN (?, ?) ORDER BY x.id LIMIT 1""",
        [plain, _canonical_dir(dir)]))
    isempty(rows) ? nothing : _named_experiment(first(rows))
end

function experiment_at_path(db, dir)
    plain = rstrip(String(dir), '/')
    rows = sqlrows(DBInterface.execute(db,
        "SELECT id, name, path FROM experiments WHERE rtrim(path, '/') IN (?, ?) ORDER BY id LIMIT 1",
        [plain, _canonical_dir(dir)]))
    isempty(rows) ? nothing : _named_experiment(first(rows))
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
         :mode => r.mode,
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

expand_folder(p::AbstractString) = (q = String(strip(p)); q == "~" || startswith(q, "~/") ? expanduser(q) : q)

function _natural_key(s::AbstractString)
    parts = Any[]
    for m in eachmatch(r"\d+|\D+", lowercase(s))
        push!(parts, isdigit(m.match[1]) ? (0, parse(BigInt, m.match), "") : (1, big(0), m.match))
    end
    parts
end

natural_sort(xs::Vector{String}) = xs[sortperm(map(_natural_key, xs))]

function describe_folder(db, dir::AbstractString)
    roots = setup_roots(dir)
    kind, setups = isempty(roots) ? (nothing, nothing) :
                   roots == [dir] ? ("setup", nothing) : ("beamtime", length(roots))
    added = kind == "beamtime" ? experiment_at_path(db, dir) :
            kind == "setup" ? something(experiment_using_dir(db, _setup_data_dir(dir)),
                                        experiment_at_path(db, dir), Some(nothing)) :
            nothing
    Dict(:path => dir, :kind => kind, :setups => setups,
         :added => added === nothing ? nothing : Dict(:id => added.id, :name => added.name))
end

function start_folder(db)
    rows = sqlrows(DBInterface.execute(db, "SELECT path FROM experiments WHERE path IS NOT NULL ORDER BY id DESC LIMIT 1"))
    last = isempty(rows) ? nothing : dirname(_trim_slash(String(first(rows).path)))
    last !== nothing && isdir(last) ? last : homedir()
end

_parent_folder(dir) = (t = _trim_slash(dir); isempty(t) ? nothing : (p = dirname(t); isempty(p) ? "/" : p))

function list_folder(db, path::AbstractString)
    dir = isempty(strip(path)) ? start_folder(db) : expand_folder(path)
    isabspath(dir) || return (400, "use a full path, starting with / or ~")
    isdir(dir) || return (400, "that folder does not exist")
    dir = _trim_slash(dir) == "" ? "/" : String(_trim_slash(dir))
    here = folder_listing(dir, _is_tif)
    in_setup = here.holds_image || ("data" in here.folders && _has_tifs(joinpath(dir, "data")))
    (200, Dict(:path => dir, :parent => _parent_folder(dir), :folders => natural_sort(here.folders),
               :inside_setup => in_setup))
end

function register_fs_routes!()
    @get "/api/fs/list" function(req::HTTP.Request)
        q = HTTP.queryparams(HTTP.URI(req.target))
        status, body = list_folder(current_db(), get(q, "path", ""))
        _json(status, status == 200 ? body : Dict(:error => body))
    end
    @get "/api/fs/describe" function(req::HTTP.Request)
        q = HTTP.queryparams(HTTP.URI(req.target))
        path = expand_folder(get(q, "path", ""))
        isabspath(path) || return _json(400, Dict(:error => "use a full path, starting with / or ~"))
        isdir(path) || return _json(400, Dict(:error => "that folder does not exist"))
        _json(200, describe_folder(current_db(), String(_trim_slash(path)) == "" ? "/" : String(_trim_slash(path))))
    end
    @get "/api/fs/preview" function(req::HTTP.Request)
        q    = HTTP.queryparams(HTTP.URI(req.target))
        path = expand_folder(get(q, "path", ""))
        isabspath(path) || return _json(400, Dict(:error => "use a full path, starting with / or ~"))
        isdir(path) || return _json(400, Dict(:error => "path is not a directory"))
        db = current_db()
        existing = experiment_at_path(db, path)
        _json(200, Dict(:path => path,
                        :existing => existing === nothing ? nothing :
                            Dict(:id => existing.id, :name => existing.name),
                        :setups => [_setup_preview(db, path, root) for root in setup_roots(path)]))
    end
    @get "/api/fs/manifest" function(req::HTTP.Request)
        q    = HTTP.queryparams(HTTP.URI(req.target))
        path = get(q, "path", "")
        isdir(path) || return _json(400, Dict(:error => "path is not a directory"))
        opt(k) = (v = get(q, k, ""); isempty(v) ? nothing : v)
        for k in ("image", "metadata", "integration")
            bad = pattern_problem(k, opt("$(k)_pattern"))
            bad === nothing || return _json(400, Dict(:error => bad))
        end
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
