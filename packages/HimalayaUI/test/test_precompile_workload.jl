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
    "POST /api/experiments", "PATCH /api/experiments/{id}", "DELETE /api/experiments/{id}",
    "POST /api/experiments/{id}/scan",
    "GET /api/experiments/{id}/export", "GET /api/experiments/{id}/samples",
    "GET /api/fs/manifest", "GET /api/fs/preview", "GET /api/fs/suggest",
    "GET /api/resolve",
    "PATCH /api/samples/{id}",
    "DELETE /api/indices/{id}",
    "POST /api/exposures/{id}/speculative",
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
