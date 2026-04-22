using HTTP, JSON3, DBInterface, Tables, Oxygen, SQLite

"""
    mark_all_indices_stale!(db, exposure_id) -> Int

Set `status = 'stale'` on every index for this exposure. Returns the count
that was marked.
"""
function mark_all_indices_stale!(db::SQLite.DB, exposure_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM indices WHERE exposure_id = ?", [exposure_id]))
    n = Int(rows[1].c)
    DBInterface.execute(db,
        "UPDATE indices SET status = 'stale' WHERE exposure_id = ?", [exposure_id])
    n
end

function register_peaks_routes!()
    @get "/api/exposures/{id}/peaks" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, exposure_id, q, intensity, prominence, sharpness, source
             FROM peaks WHERE exposure_id = ? ORDER BY q", [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(rows_to_json(rows)))
    end

    @post "/api/exposures/{id}/peaks" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        q    = Float64(body.q)

        res = DBInterface.execute(db,
            "INSERT INTO peaks (exposure_id, q, source) VALUES (?, ?, 'manual')",
            [id, q])
        peak_id = Int(DBInterface.lastrowid(res))

        stale = mark_all_indices_stale!(db, id)

        log_action!(db, req; action = "add_peak",
            entity_type = "peak", entity_id = peak_id, note = "q=$q")

        HTTP.Response(201, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => peak_id, :exposure_id => id,
                             :q => q, :source => "manual",
                             :stale_indices => stale)))
    end

    @delete "/api/peaks/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id FROM peaks WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "peak not found")))
        exposure_id = Int(rows[1].exposure_id)

        DBInterface.execute(db,
            "DELETE FROM index_peaks WHERE peak_id = ?", [id])
        DBInterface.execute(db, "DELETE FROM peaks WHERE id = ?", [id])

        mark_all_indices_stale!(db, exposure_id)

        log_action!(db, req; action = "remove_peak",
            entity_type = "peak", entity_id = id)

        HTTP.Response(204)
    end
end
