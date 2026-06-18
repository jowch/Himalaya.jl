using HTTP, DBInterface, Tables, SQLite

"""
    get_username(req) -> Union{String, Nothing}

Return the `X-Username` header value if present and non-empty, else nothing.
"""
function get_username(req::HTTP.Request)
    v = HTTP.header(req, "X-Username", "")
    isempty(v) ? nothing : String(v)
end

"""
    get_client_id(req) -> Union{String, Nothing}

Return the `X-Client-Id` header value if present and non-empty, else nothing.
This is the per-tab SSE routing identity, distinct from `X-Username` (audit
identity). See docs/event-log.md.
"""
function get_client_id(req::HTTP.Request)
    v = HTTP.header(req, "X-Client-Id", "")
    isempty(v) ? nothing : String(v)
end

"""
    get_client_op_id(req) -> Union{String, Nothing}

Return the `X-Client-Op-Id` header value if present and non-empty, else nothing.
This is the per-mutation idempotency key sent by the client (one UUID per
queued op), distinct from `X-Client-Id` (per-tab SSE routing identity).
See docs/mutation-queue.md.
"""
function get_client_op_id(req::HTTP.Request)
    v = HTTP.header(req, "X-Client-Op-Id", "")
    isempty(v) ? nothing : String(v)
end

"""
    get_or_create_user!(db, username) -> user_id::Int
"""
function get_or_create_user!(db::SQLite.DB, username::String)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM users WHERE username = ?", [username]))
    isempty(rows) || return Int(rows[1].id)
    res = DBInterface.execute(db,
        "INSERT INTO users (username) VALUES (?)", [username])
    Int(DBInterface.lastrowid(res))
end

"""
    get_user_id_for_request(db, req) -> Union{Int, Nothing}

Resolve the request's `X-Username` header to a `users.id`, creating the row
if needed (mirrors `get_or_create_user!` in the action-logging path). Returns
`nothing` when the header is absent or empty — used by route author gates
that 403 unauthenticated callers without minting a phantom user row.
"""
function get_user_id_for_request(db::SQLite.DB, req::HTTP.Request)::Union{Int, Nothing}
    username = get_username(req)
    username === nothing && return nothing
    get_or_create_user!(db, username)
end

"""
    log_action!(db, req; action, entity_type, entity_id, note=nothing)

Record an entry in `user_actions`. Missing `X-Username` => user_id = NULL.

Thin wrapper around `apply_event!` for backwards compatibility. Existing
callers continue to work unchanged. New code in R4.2 should call
`apply_event!` directly with structured payloads.
"""
function log_action!(db::SQLite.DB, req::HTTP.Request;
        action::String,
        entity_type::String,
        entity_id::Integer,
        note::Union{String, Nothing} = nothing)
    apply_event!(db, req;
                 kind        = action,
                 entity_type = entity_type,
                 entity_id   = entity_id,
                 payload     = note === nothing ? nothing : Dict(:note => note))
    nothing
end
