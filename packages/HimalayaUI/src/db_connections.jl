using SQLite, DBInterface, HTTP

mutable struct ReaderPool
    path::String
    capacity::Int
    idle::Channel{SQLite.DB}
    opened::Vector{SQLite.DB}
    lock::ReentrantLock
end

mutable struct ReaderScope
    pool::ReaderPool
    reader::Union{SQLite.DB, Nothing}
end

const _READER_POOL = Ref{Union{ReaderPool, Nothing}}(nothing)
const _READER_SCOPE_KEY = :himalaya_reader_scope
const _UNPOOLED_GET_PATHS = ("/api/events",)

_poolable(path::AbstractString) = !isempty(path) && path != ":memory:" && !startswith(path, "file:")

function _open_reader(path::AbstractString)
    reader = SQLite.DB(path)
    DBInterface.execute(reader, "PRAGMA busy_timeout = 5000")
    DBInterface.execute(reader, "PRAGMA query_only = ON")
    reader
end

_new_reader_pool(db::Union{SQLite.DB, Nothing}) =
    db !== nothing && _poolable(db.file) ?
        ReaderPool(db.file, max(4, 2 * Threads.nthreads()), Channel{SQLite.DB}(Inf), SQLite.DB[], ReentrantLock()) :
        nothing

function _close_reader_pool!(pool::Union{ReaderPool, Nothing})
    pool === nothing && return nothing
    lock(pool.lock) do
        while isready(pool.idle)
            take!(pool.idle)
        end
        foreach(r -> isopen(r) && DBInterface.close!(r), pool.opened)
        empty!(pool.opened)
    end
    nothing
end

function _reset_reader_pool!(db::Union{SQLite.DB, Nothing})
    old = _READER_POOL[]
    _READER_POOL[] = _new_reader_pool(db)
    _close_reader_pool!(old)
end

function _checkout_reader(pool::ReaderPool)
    isready(pool.idle) && return take!(pool.idle)
    fresh = lock(pool.lock) do
        length(pool.opened) < pool.capacity || return nothing
        reader = _open_reader(pool.path)
        push!(pool.opened, reader)
        reader
    end
    fresh === nothing ? take!(pool.idle) : fresh
end

function _return_reader(pool::ReaderPool, reader::SQLite.DB)
    if pool === _READER_POOL[] && isopen(reader)
        put!(pool.idle, reader)
    elseif isopen(reader)
        DBInterface.close!(reader)
    end
    nothing
end

function _scoped_reader()
    scope = get(task_local_storage(), _READER_SCOPE_KEY, nothing)
    scope === nothing && return nothing
    if scope.reader === nothing
        reader = _checkout_reader(scope.pool)
        try
            DBInterface.execute(reader, "BEGIN")
        catch
            _discard_reader!(scope.pool, reader)
            rethrow()
        end
        scope.reader = reader
    end
    scope.reader
end

function _discard_reader!(pool::ReaderPool, reader::SQLite.DB)
    lock(() -> filter!(r -> r !== reader, pool.opened), pool.lock)
    isopen(reader) && DBInterface.close!(reader)
    nothing
end

function _end_snapshot!(pool::ReaderPool, reader::SQLite.DB)
    try
        SQLite.finalize_statements!(reader)
        DBInterface.execute(reader, "ROLLBACK")
        true
    catch
        _discard_reader!(pool, reader)
        false
    end
end

function _with_reader(f)
    pool = _READER_POOL[]
    pool === nothing && return f()
    scope = ReaderScope(pool, nothing)
    try
        task_local_storage(f, _READER_SCOPE_KEY, scope)
    finally
        reader = scope.reader
        reader === nothing || (_end_snapshot!(pool, reader) && _return_reader(pool, reader))
    end
end

_is_read(req::HTTP.Request) = req.method == "GET" || req.method == "HEAD"
_is_event_stream(req::HTTP.Request) = HTTP.URI(req.target).path in _UNPOOLED_GET_PATHS

with_background_db(f, db::SQLite.DB) = _poolable(db.file) ? with_background_db(f, String(db.file)) : f(db)

function with_background_db(f, path::String)
    background = SQLite.DB(path)
    try
        DBInterface.execute(background, "PRAGMA foreign_keys = ON")
        DBInterface.execute(background, "PRAGMA busy_timeout = 10000")
        f(background)
    finally
        DBInterface.close!(background)
    end
end

function db_connection_middleware(handler)
    function (req)
        req isa HTTP.Request || return handler(req)
        _is_read(req) || return lock(() -> handler(req), _DB_WRITE_LOCK)
        _is_event_stream(req) && return handler(req)
        _with_reader(() -> handler(req))
    end
end
