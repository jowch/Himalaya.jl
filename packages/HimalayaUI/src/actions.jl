using HTTP, DBInterface, Tables, SQLite

function get_username(req::HTTP.Request)
    v = HTTP.header(req, "X-Username", "")
    isempty(v) ? nothing : String(v)
end

function get_client_id(req::HTTP.Request)
    v = HTTP.header(req, "X-Client-Id", "")
    isempty(v) ? nothing : String(v)
end

function get_client_op_id(req::HTTP.Request)
    v = HTTP.header(req, "X-Client-Op-Id", "")
    isempty(v) ? nothing : String(v)
end

function get_or_create_user!(db::SQLite.DB, username::String)
    rows = sqlrows(DBInterface.execute(db,
        "SELECT id FROM users WHERE username = ?", [username]))
    isempty(rows) || return Int(rows[1].id)
    res = DBInterface.execute(db,
        "INSERT INTO users (username) VALUES (?)", [username])
    Int(DBInterface.lastrowid(res))
end

function get_user_id_for_request(db::SQLite.DB, req::HTTP.Request)::Union{Int, Nothing}
    username = get_username(req)
    username === nothing && return nothing
    get_or_create_user!(db, username)
end

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
