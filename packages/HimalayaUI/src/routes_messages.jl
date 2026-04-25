using HTTP, JSON3, DBInterface, Tables, Oxygen, SQLite

"""
    select_message_sql(where_clause)

Build a `SELECT … LEFT JOIN users …` query that returns the canonical message
row shape used by both GET (list) and POST (insert echo):
`{ id, sample_id, author_id, author, body, created_at }`.
`author` is the current username at read-time, or `null` if the user was
deleted (`author_id` is preserved as the FK).
"""
function select_message_sql(where_clause::String)
    """
    SELECT m.id, m.sample_id, m.author_id,
           u.username AS author,
           m.body, m.created_at
    FROM sample_messages m
    LEFT JOIN users u ON u.id = m.author_id
    $where_clause
    """
end

function register_messages_routes!()
    @get "/api/samples/{id}/messages" function(req::HTTP.Request, id::Int)
        db   = current_db()
        sql  = select_message_sql(
            "WHERE m.sample_id = ? ORDER BY m.created_at ASC, m.id ASC")
        rows = Tables.rowtable(DBInterface.execute(db, sql, [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(rows_to_json(rows)))
    end

    @post "/api/samples/{id}/messages" function(req::HTTP.Request, id::Int)
        db       = current_db()
        username = get_username(req)
        if username === nothing
            return HTTP.Response(401,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "X-Username header required")))
        end

        body = json(req)
        text = haskey(body, :body) ? strip(String(body.body)) : ""
        if isempty(text)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "message body required")))
        end

        author_id = get_or_create_user!(db, username)
        res = DBInterface.execute(db,
            "INSERT INTO sample_messages (sample_id, author_id, body) VALUES (?, ?, ?)",
            [id, author_id, text])
        msg_id = Int(DBInterface.lastrowid(res))

        sql = select_message_sql("WHERE m.id = ?")
        row = Tables.rowtable(DBInterface.execute(db, sql, [msg_id]))[1]

        log_action!(db, req; action = "add_message",
            entity_type = "sample_message", entity_id = msg_id)

        HTTP.Response(201, ["Content-Type" => "application/json"],
            JSON3.write(row_to_json(row)))
    end
end
