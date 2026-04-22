using HTTP, JSON3, DBInterface, Tables, Oxygen, SQLite

"""
    ensure_custom_group!(db, exposure_id) -> (group_id, created)

Returns the id of the custom group for this exposure and whether it was just
created. If a custom group already exists, returns (existing_id, false).
Otherwise clones the auto group's members into a new custom group (active),
demotes the auto group to inactive, returns (new_id, true).

Errors if no auto group exists for the exposure.
"""
function ensure_custom_group!(db::SQLite.DB, exposure_id::Int)
    existing = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM index_groups
         WHERE exposure_id = ? AND kind = 'custom'", [exposure_id]))
    isempty(existing) || return (Int(existing[1].id), false)

    auto_rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM index_groups
         WHERE exposure_id = ? AND kind = 'auto'", [exposure_id]))
    isempty(auto_rows) && error("no auto group for exposure $exposure_id")
    auto_id = Int(auto_rows[1].id)

    res = DBInterface.execute(db,
        "INSERT INTO index_groups (exposure_id, kind, active)
         VALUES (?, 'custom', 1)", [exposure_id])
    custom_id = Int(DBInterface.lastrowid(res))

    DBInterface.execute(db,
        "INSERT INTO index_group_members (group_id, index_id)
         SELECT ?, index_id FROM index_group_members WHERE group_id = ?",
        [custom_id, auto_id])

    DBInterface.execute(db,
        "UPDATE index_groups SET active = 0 WHERE id = ?", [auto_id])

    (custom_id, true)
end

function _group_with_members(db::SQLite.DB, group_id::Int)
    g = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM index_groups WHERE id = ?", [group_id]))[1]
    members = Tables.rowtable(DBInterface.execute(db,
        "SELECT index_id FROM index_group_members
         WHERE group_id = ? ORDER BY index_id", [group_id]))
    d = row_to_json(g; bool_keys = (:active,))
    d[:members] = [Int(m.index_id) for m in members]
    d
end

function register_analysis_routes!()
    @get "/api/exposures/{id}/indices" function(req::HTTP.Request, id::Int)
        db = current_db()
        indices = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM indices WHERE exposure_id = ? ORDER BY score DESC", [id]))
        out = map(indices) do ix
            peaks = Tables.rowtable(DBInterface.execute(db,
                "SELECT peak_id, ratio_position, residual
                 FROM index_peaks WHERE index_id = ? ORDER BY ratio_position",
                [Int(ix.id)]))
            d = row_to_json(ix)
            d[:peaks] = rows_to_json(peaks)
            d
        end
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(out))
    end

    @get "/api/exposures/{id}/groups" function(req::HTTP.Request, id::Int)
        db     = current_db()
        groups = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM index_groups WHERE exposure_id = ? ORDER BY id", [id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write([_group_with_members(db, Int(g.id)) for g in groups]))
    end

    @post "/api/groups/{id}/members" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        index_id = Int(body.index_id)

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id, kind FROM index_groups WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "group not found")))
        exposure_id = Int(rows[1].exposure_id)

        custom_id, _ = ensure_custom_group!(db, exposure_id)

        DBInterface.execute(db,
            "INSERT OR IGNORE INTO index_group_members (group_id, index_id)
             VALUES (?, ?)", [custom_id, index_id])

        log_action!(db, req; action = "confirm_index",
            entity_type = "index", entity_id = index_id)

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(_group_with_members(db, custom_id)))
    end

    @delete "/api/groups/{id}/members/{index_id}" function(req::HTTP.Request, id::Int, index_id::Int)
        db = current_db()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT exposure_id FROM index_groups WHERE id = ?", [id]))
        isempty(rows) && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "group not found")))
        exposure_id = Int(rows[1].exposure_id)

        custom_id, _ = ensure_custom_group!(db, exposure_id)
        DBInterface.execute(db,
            "DELETE FROM index_group_members
             WHERE group_id = ? AND index_id = ?", [custom_id, index_id])

        log_action!(db, req; action = "exclude_index",
            entity_type = "index", entity_id = index_id)

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(_group_with_members(db, custom_id)))
    end
end
