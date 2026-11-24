const _STATE_CONTAINER = Union{AbstractDict, AbstractVector, AbstractSet}

function _module_state()
    state = Dict{Symbol, Any}()
    for name in names(HimalayaUI; all = true)
        isdefined(HimalayaUI, name) && isconst(HimalayaUI, name) || continue
        value = getfield(HimalayaUI, name)
        if value isa _STATE_CONTAINER
            state[name] = copy(value)
        elseif value isa Base.RefValue
            held = value[]
            state[name] = held isa _STATE_CONTAINER ? copy(held) : held
        end
    end
    state
end

const _WORKLOAD_REQUESTS = Tuple{String, String}[]

@testset "precompile workload drives every hot route and leaves no module state behind" begin
    HimalayaUI._clear_peak_set_hash_cache!()
    before = _module_state()
    @test all(haskey(before, n) for n in (:_DB_REF, :SSE_SUBSCRIBERS, :OP_LOCKS,
                                          :_PEAK_SET_HASH_CACHE, :_PEAK_SET_HASH_CACHE_ORDER))
    mktempdir() do tmp
        @test HimalayaUI.run_route_workload(tmp; record = _WORKLOAD_REQUESTS) === nothing
    end
    after = _module_state()
    @test sort!(collect(keys(after))) == sort!(collect(keys(before)))
    @test [name for name in keys(before) if !isequal(before[name], after[name])] == Symbol[]
    probe = Oxygen.internalrequest(HTTP.Request("GET", "/api/health"); serialize = true, catch_errors = true)
    @test probe.status == 404
end

const _COLD_ROUTES_OFF_THE_TRIAGE_PATH = Set([
    "GET /**",
    "POST /api/experiments", "DELETE /api/experiments/{id}",
    "POST /api/experiments/{id}/scan", "POST /api/experiments/{id}/setups",
    "GET /api/fs/manifest", "GET /api/fs/preview",
    "PATCH /api/samples/{id}",
])

function _registered_routes()
    src = joinpath(pkgdir(HimalayaUI), "src")
    routes = Set{String}()
    for file in readdir(src; join = true)
        endswith(file, ".jl") || continue
        for m in eachmatch(r"@(get|post|patch|put|delete) \"([^\"]+)\"", read(file, String))
            push!(routes, uppercase(m.captures[1]) * " " * m.captures[2])
        end
    end
    routes
end

function _route_matches(route::AbstractString, method::AbstractString, target::AbstractString)
    verb, template = split(route, ' '; limit = 2)
    pattern = Regex("^" * replace(template, "*" => "\\*", r"\{[^}]+\}" => "[^/]+") * "\$")
    verb == method && occursin(pattern, first(split(target, '?')))
end

@testset "every registered route is either in the precompile workload or listed as off the triage path" begin
    @test !isempty(_WORKLOAD_REQUESTS)
    routes = _registered_routes()
    @test length(routes) > 40
    driven = Set(r for r in routes if any(((m, t),) -> _route_matches(r, m, t), _WORKLOAD_REQUESTS))
    @test sort!(collect(setdiff(routes, driven, _COLD_ROUTES_OFF_THE_TRIAGE_PATH))) == String[]
    @test sort!(collect(intersect(driven, _COLD_ROUTES_OFF_THE_TRIAGE_PATH))) == String[]
    @test sort!(collect(setdiff(_COLD_ROUTES_OFF_THE_TRIAGE_PATH, routes))) == String[]
end

@testset "the workload's mixed trace indexes every phase the pipeline searches" begin
    mktempdir() do dir
        analysis_dir = joinpath(dir, "analysis")
        mkpath(analysis_dir)
        HimalayaUI._workload_trace(
            joinpath(analysis_dir, replace(HimalayaUI.DEFAULT_INTEGRATION_PATTERN, "{name}" => "frame_mixed")),
            HimalayaUI._workload_mixed_centers())
        db = HimalayaUI.open_db(joinpath(dir, "himalaya.db"))
        try
            experiment = HimalayaUI.init_experiment!(db; path = dir, data_dir = dir, analysis_dir)
            sample = HimalayaUI.create_sample!(db; experiment_id = experiment, name = "mixed")
            exposure = HimalayaUI.create_exposure!(db; experiment_id = experiment, sample_id = sample,
                                                   filename = "frame_mixed", image_path = joinpath(dir, "frame.tif"))
            HimalayaUI.analyze_exposure!(db, exposure, analysis_dir)
            found = Set(String(r.phase) for r in HimalayaUI.sqlrows(DBInterface.execute(db,
                "SELECT phase FROM indices WHERE exposure_id = ?", [exposure])))
            searched = Set(string(nameof(P)) for P in HimalayaUI._WORKLOAD_PHASES if P !== HimalayaUI.Himalaya.Fm3m)
            @test setdiff(searched, found) == Set{String}()
        finally
            DBInterface.close!(db)
        end
    end
end

@testset "DB reads go through sqlrows so their row types do not depend on the data" begin
    src = joinpath(pkgdir(HimalayaUI), "src")
    offenders = String[]
    for file in readdir(src; join = true)
        endswith(file, ".jl") || continue
        for (n, line) in enumerate(eachline(file))
            occursin(r"\browtable\b", line) && push!(offenders, "$(basename(file)):$n")
        end
    end
    @test offenders == String[]
    isempty(offenders) || @info """
        Read query results with `sqlrows(DBInterface.execute(...))`, not `Tables.rowtable`. \
        `Tables.rowtable` types each column by the values it finds, so a NULL or an empty result \
        makes a new row type, and the code that reads it compiles again on first use.""" offenders
end
