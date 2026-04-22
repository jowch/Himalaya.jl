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
    log_action!(db, req; action, entity_type, entity_id, note=nothing)

Record an entry in `user_actions`. Missing `X-Username` => user_id = NULL.
"""
function log_action!(db::SQLite.DB, req::HTTP.Request;
        action::String,
        entity_type::String,
        entity_id::Integer,
        note::Union{String, Nothing} = nothing)
    username = get_username(req)
    user_id  = username === nothing ? nothing : get_or_create_user!(db, username)
    DBInterface.execute(db,
        "INSERT INTO user_actions (user_id, action, entity_type, entity_id, note)
         VALUES (?, ?, ?, ?, ?)",
        [user_id, action, entity_type, Int(entity_id), note])
    nothing
end
