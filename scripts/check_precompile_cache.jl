using HimalayaUI

function _kept(mi::Core.MethodInstance)
    ci = isdefined(mi, :cache) ? mi.cache : nothing
    while ci !== nothing
        ci.max_world == typemax(UInt) && return true
        ci = isdefined(ci, :next) ? ci.next : nothing
    end
    false
end

function check(f)
    total = 0
    rejected = String[]
    for m in methods(f), mi in Base.specializations(m)
        mi isa Core.MethodInstance || continue
        total += 1
        _kept(mi) || push!(rejected, first(sprint(show, mi.specTypes), 300))
    end
    total, rejected
end

total, rejected = check(HimalayaUI.with_idempotency)
if total == 0
    println(stderr, """
        No precompiled `with_idempotency` code to check, so the route warm-up did not run.
        Delete packages/HimalayaUI/LocalPreferences.toml (or set precompile_workload = true) and precompile again.""")
    exit(1)
end
if !isempty(rejected)
    println(stderr, """
        $(length(rejected)) of $total precompiled `with_idempotency` specializations were thrown away when \
        HimalayaUI loaded, so those edit routes compile again on their first use.
        A method added after HimalayaUI was compiled invalidated something they call. The usual cause is a \
        call on an abstract argument type, such as `JSON3.read(::AbstractString)`: pass a concrete type.
        Find the culprit with SnoopCompile: `invalidation_trees(@snoop_invalidations using HimalayaUI)`.
        Rejected:""")
    foreach(r -> println(stderr, "  ", r), rejected)
    exit(1)
end
println("All $total precompiled with_idempotency specializations survived loading.")
