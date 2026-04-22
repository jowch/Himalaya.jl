using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_users_routes!()
    @get "/api/users" function(req::HTTP.Request)
        rows = Tables.rowtable(DBInterface.execute(current_db(),
            "SELECT id, username FROM users ORDER BY id"))
        rows_to_json(rows)
    end

    @post "/api/users" function(req::HTTP.Request)
        body     = json(req)
        username = String(body.username)
        isempty(username) && return HTTP.Response(400,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "username required")))

        db    = current_db()
        rows  = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, username FROM users WHERE username = ?", [username]))
        if !isempty(rows)
            return HTTP.Response(200, ["Content-Type" => "application/json"],
                JSON3.write(row_to_json(rows[1])))
        end

        uid = get_or_create_user!(db, username)
        HTTP.Response(201, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => uid, :username => username)))
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
