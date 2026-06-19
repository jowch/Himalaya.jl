# In-process route dispatch fixture — the fast analogue of with_test_server.
# Dispatches HTTP.Requests straight through Oxygen.internalrequest (no socket),
# against the same global Oxygen CONTEXT[] that register_routes! populates.
using HTTP, SQLite
import Oxygen
using HimalayaUI

# Router-liveness check (spec blocker B1). resetstate() — called by BOTH
# start_test_server! and stop_test_server! — swaps CONTEXT[] for a fresh
# empty-router ServerContext, so a wire keeper running before an in-process
# file under GROUP=All wipes the router. Probe a route register_routes! always
# mounts (/api/health) and re-register only when it's gone. Never a sticky bool.
function _ensure_inproc_routes!()
    probe = Oxygen.internalrequest(HTTP.Request("GET", "/api/health");
                                   serialize = true, catch_errors = true)
    probe.status == 200 || HimalayaUI.register_routes!()
    nothing
end

"""
    with_inproc_routes(f, db)

Bind `db` as the singleton, ensure routes are registered against the live
Oxygen context, and pass `f` a `call(method, target; headers, body)` closure
that returns the `HTTP.Response` from `Oxygen.internalrequest` (serialize=true
so thrown handlers become 4xx/5xx exactly as production serve() does).

Mirrors `stop_test_server!`'s scrub in `finally` (minus resetstate()): nulls
`_DB_REF`, clears `SSE_SUBSCRIBERS` under `SSE_LOCK`, empties `OP_LOCKS`.
"""
function with_inproc_routes(f, db::SQLite.DB)
    _ensure_inproc_routes!()
    HimalayaUI.bind_db!(db)
    call = function (method::AbstractString, target::AbstractString;
                     headers = Pair{String,String}[], body = UInt8[])
        req = HTTP.Request(method, target, headers, body)
        return Oxygen.internalrequest(req; serialize = true, catch_errors = true)
    end
    try
        f(call)
    finally
        HimalayaUI._DB_REF[] = nothing
        lock(HimalayaUI.SSE_LOCK) do
            HimalayaUI.SSE_SUBSCRIBERS[] = []
        end
        empty!(HimalayaUI.OP_LOCKS)
    end
end
