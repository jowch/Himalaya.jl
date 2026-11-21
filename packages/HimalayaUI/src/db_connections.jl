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
        scope.reader = _checkout_reader(scope.pool)
    end
    scope.reader
end

function _with_reader(f)
    pool = _READER_POOL[]
    pool === nothing && return f()
    scope = ReaderScope(pool, nothing)
    try
        task_local_storage(f, _READER_SCOPE_KEY, scope)
    finally
        scope.reader === nothing || _return_reader(pool, scope.reader)
    end
end

function _reads_through_pool(req::HTTP.Request)
    req.method == "GET" || req.method == "HEAD" || return false
    path = HTTP.URI(req.target).path
    !(path in _UNPOOLED_GET_PATHS)
end

"""
    with_background_db(f, db)

Run `f` with a connection of its own to `db`'s file, for work spawned off a
request or a timer. Such work outlives the request that started it, so it can
not share the singleton, which only the request holding `_DB_WRITE_LOCK` may
touch. Its transactions still take `_DB_WRITE_LOCK`. An in-memory `db` has no
file to reopen, so `f` gets `db` itself.
"""
function with_background_db(f, db::SQLite.DB)
    _poolable(db.file) || return f(db)
    background = SQLite.DB(db.file)
    try
        DBInterface.execute(background, "PRAGMA foreign_keys = ON")
        DBInterface.execute(background, "PRAGMA busy_timeout = 10000")
        f(background)
    finally
        DBInterface.close!(background)
    end
end

"""
    db_connection_middleware(handler)

Oxygen middleware that gives each request a connection that is safe for its
kind. GET and HEAD requests read through a pooled, `query_only` reader
connection, so they never touch the singleton's statement cache and never wait
for a writer (WAL). Every other request holds `_DB_WRITE_LOCK` for its whole
duration, so the singleton is only ever used by one task at a time. The event
stream keeps the singleton because it holds its request open indefinitely and
does not read the DB.
"""
function db_connection_middleware(handler)
    function (req)
        req isa HTTP.Request || return handler(req)
        _reads_through_pool(req) && return _with_reader(() -> handler(req))
        lock(() -> handler(req), _DB_WRITE_LOCK)
    end
end
