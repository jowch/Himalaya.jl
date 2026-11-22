using SQLite, DBInterface, HTTP, Tables

const OP_LOCKS    = Dict{String, ReentrantLock}()
const OP_LOCKS_MU = ReentrantLock()

function _op_lock(op_id::String)::ReentrantLock
    lock(OP_LOCKS_MU) do
        get!(OP_LOCKS, op_id, ReentrantLock())
    end
end

function _drop_op_lock!(op_id::String)
    lock(OP_LOCKS_MU) do
        delete!(OP_LOCKS, op_id)
    end
end

function _lookup_cached_response(db::SQLite.DB, op_id::String)::Union{HTTP.Response, Nothing}
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT status_code, body FROM idempotent_responses WHERE client_op_id = ?",
        [op_id]))
    isempty(rows) && return nothing
    row = rows[1]
    return HTTP.Response(Int(row.status_code),
                         ["Content-Type" => "application/json"];
                         body = String(row.body))
end

function _write_transaction(f, db::SQLite.DB)
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(f, db)
    end
end

function with_idempotency(f, db::SQLite.DB, req::HTTP.Request)
    op_id = get_client_op_id(req)
    if op_id === nothing
        try
            response = _write_transaction(f, db)
            _flush_post_commit_broadcasts!()
            return response
        catch
            _clear_post_commit_broadcasts!()
            rethrow()
        end
    end

    cached = _lookup_cached_response(db, op_id)
    cached === nothing || return cached

    return lock(_op_lock(op_id)) do
        replayed_cache = false
        local response
        try
            response = _write_transaction(db) do
                cached2 = _lookup_cached_response(db, op_id)
                if cached2 !== nothing
                    replayed_cache = true
                    return cached2
                end
                resp = f()
                if resp.status < 400
                    DBInterface.execute(db,
                        "INSERT INTO idempotent_responses (client_op_id, status_code, body) VALUES (?, ?, ?)",
                        [op_id, Int(resp.status), String(copy(resp.body))])
                end
                return resp
            end
        catch
            _clear_post_commit_broadcasts!()
            _drop_op_lock!(op_id)
            rethrow()
        end

        if replayed_cache
            _clear_post_commit_broadcasts!()
        elseif response.status >= 400
            _clear_post_commit_broadcasts!()
            _drop_op_lock!(op_id)
        else
            _flush_post_commit_broadcasts!()
        end
        return response
    end
end

function gc_idempotent_responses!(db::SQLite.DB; ttl_seconds::Int = 3600)
    expired = lock(_DB_WRITE_LOCK) do
        Tables.rowtable(DBInterface.execute(db,
            """DELETE FROM idempotent_responses
               WHERE created_at < datetime('now', ?)
               RETURNING client_op_id""",
            ["-$(ttl_seconds) seconds"]))
    end
    isempty(expired) && return nothing
    lock(OP_LOCKS_MU) do
        for r in expired
            delete!(OP_LOCKS, String(r.client_op_id))
        end
    end
    nothing
end
