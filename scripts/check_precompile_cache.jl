using HimalayaUI
using Oxygen

function _kept(mi::Core.MethodInstance)
    ci = isdefined(mi, :cache) ? mi.cache : nothing
    while ci !== nothing
        ci.max_world == typemax(UInt) && return true
        ci = isdefined(ci, :next) ? ci.next : nothing
    end
    false
end

_checked(m::Method) =
    (m.module === HimalayaUI && !endswith(String(m.file), "precompile.jl")) ||
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
        $(length(rejected)) of $total precompiled HimalayaUI specializations were thrown away when HimalayaUI \
        loaded, so the routes that use them compile again on first use.
        A method added after HimalayaUI was compiled invalidated something they call. The usual cause is a \
        call on an abstract argument type, such as `JSON3.read(::AbstractString)` or `replace(::AbstractString, …)`: \
        convert to `String` first. A route that returns a plain Julia value goes through Oxygen's JSON.jl \
        serializer, which is invalidated too: return an `HTTP.Response` with `JSON3.write`.
        Find the culprit with SnoopCompile: `invalidation_trees(@snoop_invalidations using HimalayaUI)`.
        Rejected:""")
    foreach(r -> println(stderr, "  ", r), rejected)
    exit(1)
end
println("All $total precompiled HimalayaUI specializations survived loading.")
