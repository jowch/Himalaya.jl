using HTTP, DBInterface, Tables, Oxygen

function register_users_routes!()
    @get "/api/users" function(req::HTTP.Request)
        rows = sqlrows(DBInterface.execute(current_db(),
            "SELECT id, username, first_name, last_name FROM users ORDER BY id"))
        _json(200, rows_to_json(rows))
    end

    @post "/api/users" function(req::HTTP.Request)
        body = json(req)
        if !haskey(body, :username)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                to_json(Dict(:error => "missing field: username")))
        end
        username = String(body.username)
        isempty(username) && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            to_json(Dict(:error => "username required")))
        occursin(r"^[A-Za-z0-9_]+$", username) || return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            to_json(Dict(:error => "username must contain only letters, digits, and underscores")))

        first_name = haskey(body, :first_name) ? (isnothing(body.first_name) ? nothing : String(body.first_name)) : nothing
        last_name  = haskey(body, :last_name)  ? (isnothing(body.last_name)  ? nothing : String(body.last_name))  : nothing

        db   = current_db()
        rows = sqlrows(DBInterface.execute(db,
            "SELECT id, username, first_name, last_name FROM users WHERE username = ?", [username]))
        if !isempty(rows)
            existing = rows[1]
            existing_first = ismissing(existing.first_name) ? nothing : existing.first_name
            existing_last  = ismissing(existing.last_name)  ? nothing : existing.last_name
            new_first = (existing_first === nothing && first_name !== nothing) ? first_name : existing_first
            new_last  = (existing_last  === nothing && last_name  !== nothing) ? last_name  : existing_last
            if new_first !== existing_first || new_last !== existing_last
                DBInterface.execute(db,
                    "UPDATE users SET first_name = ?, last_name = ? WHERE id = ?",
                    [new_first, new_last, existing.id])
                rows = sqlrows(DBInterface.execute(db,
                    "SELECT id, username, first_name, last_name FROM users WHERE id = ?", [existing.id]))
            end
            return HTTP.Response(200, ["Content-Type" => "application/json"],
                to_json(row_to_json(rows[1])))
        end

        result = DBInterface.execute(db,
            "INSERT INTO users (username, first_name, last_name) VALUES (?, ?, ?)",
            [username, first_name, last_name])
        uid = Int(DBInterface.lastrowid(result))
        HTTP.Response(201, ["Content-Type" => "application/json"],
            to_json(Dict(:id => uid, :username => username,
                             :first_name => first_name, :last_name => last_name)))
    end
end
