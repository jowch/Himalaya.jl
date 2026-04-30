using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_users_routes!()
    @get "/api/users" function(req::HTTP.Request)
        rows = Tables.rowtable(DBInterface.execute(current_db(),
            "SELECT id, username, first_name, last_name FROM users ORDER BY id"))
        rows_to_json(rows)
    end

    @post "/api/users" function(req::HTTP.Request)
        body     = json(req)
        username = String(body.username)
        isempty(username) && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "username required")))
        occursin(r"^[A-Za-z0-9_]+$", username) || return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "username must contain only letters, digits, and underscores")))

        first_name = haskey(body, :first_name) ? (isnothing(body.first_name) ? nothing : String(body.first_name)) : nothing
        last_name  = haskey(body, :last_name)  ? (isnothing(body.last_name)  ? nothing : String(body.last_name))  : nothing

        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, username, first_name, last_name FROM users WHERE username = ?", [username]))
        if !isempty(rows)
            existing = rows[1]
            # Idempotent enrichment: when a previous create stored NULL for first/last,
            # a follow-up POST that supplies them should fill the gap. Never overwrite
            # an existing non-null name — preserves "first wins" for actual identity.
            existing_first = ismissing(existing.first_name) ? nothing : existing.first_name
            existing_last  = ismissing(existing.last_name)  ? nothing : existing.last_name
            new_first = (existing_first === nothing && first_name !== nothing) ? first_name : existing_first
            new_last  = (existing_last  === nothing && last_name  !== nothing) ? last_name  : existing_last
            if new_first !== existing_first || new_last !== existing_last
                DBInterface.execute(db,
                    "UPDATE users SET first_name = ?, last_name = ? WHERE id = ?",
                    [new_first, new_last, existing.id])
                rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, username, first_name, last_name FROM users WHERE id = ?", [existing.id]))
            end
            return HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(row_to_json(rows[1])))
        end

        result = DBInterface.execute(db,
            "INSERT INTO users (username, first_name, last_name) VALUES (?, ?, ?)",
            [username, first_name, last_name])
        uid = Int(DBInterface.lastrowid(result))
        HTTP.Response(201, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => uid, :username => username,
                             :first_name => first_name, :last_name => last_name)))
    end

    @get "/api/users/{username}/actions" function(req::HTTP.Request, username::String)
        db = current_db()
        urows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM users WHERE username = ?", [username]))
        isempty(urows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "user not found")))

        uid  = Int(urows[1].id)
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, timestamp, action, entity_type, entity_id, note
             FROM user_actions WHERE user_id = ? ORDER BY id DESC", [uid]))
        rows_to_json(rows)
    end
end
