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
            "SELECT id, exposure_id, q, intensity, prominence, sharpness, source, excluded
             FROM peaks WHERE exposure_id = ? ORDER BY q", [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(rows_to_json(rows; bool_keys = (:excluded,))))
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
                             :q => q, :source => "manual", :excluded => false,
                             :stale_indices => stale)))
    end

    @patch "/api/peaks/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, source FROM peaks WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "peak not found")))
        peak_row     = rows[1]
        exposure_id  = Int(peak_row.exposure_id)
        source       = String(peak_row.source)

        # Currently the only patchable field is `excluded`, and only on auto peaks.
        if !haskey(body, :excluded)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "missing field: excluded")))
        end
        if source != "auto"
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "only auto peaks can be excluded; delete manual peaks instead")))
        end

        excluded = Bool(body.excluded)
        DBInterface.execute(db,
            "UPDATE peaks SET excluded = ? WHERE id = ?", [Int(excluded), id])

        stale = mark_all_indices_stale!(db, exposure_id)

        log_action!(db, req; action = excluded ? "exclude_peak" : "include_peak",
            entity_type = "peak", entity_id = id)

        # Return the updated row
        out = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, exposure_id, q, intensity, prominence, sharpness, source, excluded
             FROM peaks WHERE id = ?", [id]))[1]

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(merge(row_to_json(out; bool_keys = (:excluded,)),
                              Dict(:stale_indices => stale))))
    end

    @delete "/api/peaks/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, source FROM peaks WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "peak not found")))
        exposure_id = Int(rows[1].exposure_id)
        source      = String(rows[1].source)

        # Manual peaks can be hard-deleted (the user owns them). Auto peaks should
        # use PATCH { excluded: true } instead — deleting them risks the next
        # reanalysis re-creating them and the user wondering why their override
        # vanished. We allow it for backwards compatibility but log a note.
        if source == "auto"
            log_action!(db, req; action = "delete_auto_peak_legacy",
                entity_type = "peak", entity_id = id,
                note = "consider PATCH { excluded: true } instead")
        end

        DBInterface.execute(db,
            "DELETE FROM index_peaks WHERE peak_id = ?", [id])
        DBInterface.execute(db, "DELETE FROM peaks WHERE id = ?", [id])

        mark_all_indices_stale!(db, exposure_id)

        log_action!(db, req; action = "remove_peak",
            entity_type = "peak", entity_id = id)

        HTTP.Response(204)
    end
end
