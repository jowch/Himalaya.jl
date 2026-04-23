using Oxygen
using HTTP
using SQLite
using JSON3
import Sockets

const _DB_REF = Ref{Union{SQLite.DB, Nothing}}(nothing)

function current_db()
    db = _DB_REF[]
    db === nothing && error("no DB bound; call serve(db; ...) or start_test_server!")
    db
end

function bind_db!(db::SQLite.DB)
    _DB_REF[] = db
    nothing
end

function register_routes!()
    # Static frontend assets — only mounted if dist directory exists with content
    dist_dir = joinpath(pkgdir(HimalayaUI), "frontend", "dist")
    if isdir(dist_dir)
        Oxygen.dynamicfiles(dist_dir, "/")
    end

    @get "/api/health" function(req::HTTP.Request)
        Dict("status" => "ok")
    end
    register_users_routes!()
    register_experiments_routes!()
    register_samples_routes!()
    register_exposures_routes!()
    register_peaks_routes!()
    register_trace_routes!()
    register_analysis_routes!()
    register_export_routes!()
end

function serve(db::SQLite.DB; host::String = "127.0.0.1", port::Int = 8080)
    Oxygen.resetstate()
    bind_db!(db)
    register_routes!()
    Oxygen.serve(; host, port, show_banner = false, docs = false, metrics = false)
end

function start_test_server!(db::SQLite.DB, port::Int)
    Oxygen.resetstate()
    bind_db!(db)
    register_routes!()
    Oxygen.serve(; host = "127.0.0.1", port, async = true, show_banner = false, docs = false, metrics = false)
    _wait_for_server(port)
end

function stop_test_server!()
    Oxygen.terminate()
    Oxygen.resetstate()
    _DB_REF[] = nothing
    nothing
end

function find_free_port()
    server = Sockets.listen(Sockets.IPv4(0), 0)
    port   = Sockets.getsockname(server)[2]
    close(server)
    Int(port)
end

function _wait_for_server(port::Int; timeout_s = 5.0)
    deadline = time() + timeout_s
    while time() < deadline
        try
            resp = HTTP.get("http://127.0.0.1:$port/api/health";
                            connect_timeout = 1, readtimeout = 1,
                            retry = false, status_exception = false)
            resp.status == 200 && return
        catch
        end
        sleep(0.05)
    end
    error("test server on port $port did not become ready within $(timeout_s)s")
end
