using HimalayaUI
using Oxygen

_live(ci::Core.CodeInstance) = ci.max_world == typemax(UInt)

function _cached(mi::Core.MethodInstance)
    ci = isdefined(mi, :cache) ? mi.cache : nothing
    while ci !== nothing
        _live(ci) && return true
        ci = isdefined(ci, :next) ? ci.next : nothing
    end
    false
end

_in_workload(mi) = mi isa Core.MethodInstance && mi.def isa Method && endswith(String(mi.def.file), "precompile.jl")

_live_caller(b) =
    (b isa Core.CodeInstance && !_in_workload(b.def) && _live(b)) ||
    (b isa Core.MethodInstance && !_in_workload(b) && _cached(b))

_called_from_live_code(mi::Core.MethodInstance) = isdefined(mi, :backedges) && any(_live_caller, mi.backedges)

_kept(mi::Core.MethodInstance) = _cached(mi) || _called_from_live_code(mi)

const _MODULES = (HimalayaUI, HimalayaUI.Himalaya, HimalayaUI.HimalayaPeakModel)

_checked(m::Method) =
    (parentmodule(m) in _MODULES && !endswith(String(m.file), "precompile.jl")) ||
    (parentmodule(m) === Oxygen.Core.Util && m.name === :format_response!)

function check()
    total = 0
    idempotent = 0
    rejected = String[]
    Base.visit(Core.methodtable) do m
        _checked(m) || return
        for mi in Base.specializations(m)
            mi isa Core.MethodInstance || continue
            Base.isdispatchtuple(mi.specTypes) || continue
            total += 1
            m.name === :with_idempotency && (idempotent += 1)
            _kept(mi) || push!(rejected, first(sprint(show, mi.specTypes), 300))
        end
    end
    total, idempotent, rejected
end

total, idempotent, rejected = check()
if idempotent == 0
    println(stderr, """
        No precompiled `with_idempotency` code to check, so the route warm-up did not run.
        Delete packages/HimalayaUI/LocalPreferences.toml (or set precompile_workload = true) and precompile again.""")
    exit(1)
end
if !isempty(rejected)
    println(stderr, """
        $(length(rejected)) of $total precompiled specializations of HimalayaUI, Himalaya and HimalayaPeakModel were thrown away when \
        HimalayaUI loaded, so the routes that use them compile again on first use.
        A method added after HimalayaUI was compiled invalidated something they call. The usual cause is a \
        call on an abstract argument type, such as `JSON.parse(::AbstractString)` or `replace(::AbstractString, …)`: \
        convert to `String` first. A route that returns a plain Julia value goes through Oxygen's JSON.jl \
        serializer, which is invalidated too: return an `HTTP.Response` with `to_json`, never `JSON.json`.
        Find the culprit with SnoopCompile: `invalidation_trees(@snoop_invalidations using HimalayaUI)`.
        Rejected:""")
    foreach(r -> println(stderr, "  ", r), rejected)
    exit(1)
end
println("All $total precompiled specializations of HimalayaUI, Himalaya and HimalayaPeakModel survived loading.")
