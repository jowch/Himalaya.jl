using HTTP, JSON3, DBInterface, Tables, Oxygen, SQLite

"""
    _peak_curation_tol(q) -> Float64

Tolerance for matching auto peak q-values to exclude-curation q-values.
"""
_peak_curation_tol(q::Real) = max(1e-6, abs(Float64(q)) * 0.001)

"""
    _effective_peak_row(db, id) -> Union{NamedTuple, Nothing}

Look up a peak by id across both auto_peaks and peak_curations.
Returns a NamedTuple with fields: id, exposure_id, q, intensity, prominence,
sharpness, source, excluded — or nothing if not found.
"""
function _effective_peak_row(db::SQLite.DB, id::Int)
    # Try auto_peaks first
    auto_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT a.id, a.exposure_id, a.q, a.intensity, a.prominence, a.sharpness,
                  'auto' AS source,
                  CASE WHEN c.q IS NOT NULL THEN 1 ELSE 0 END AS excluded
           FROM auto_peaks a
           LEFT JOIN peak_curations c
               ON c.exposure_id = a.exposure_id
              AND c.kind = 'exclude'
              AND ABS(c.q - a.q) <= MAX(1e-6, ABS(a.q) * 0.001)
           WHERE a.id = ?""", [id]))
    isempty(auto_rows) || return auto_rows[1]

    # Try peak_curations (kind='add')
    curation_rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT id, exposure_id, q, NULL AS intensity, NULL AS prominence,
                  NULL AS sharpness, 'manual' AS source, 0 AS excluded
           FROM peak_curations WHERE id = ? AND kind = 'add'""", [id]))
    isempty(curation_rows) && return nothing
    return curation_rows[1]
end

function register_peaks_routes!()
    @get "/api/exposures/{id}/peaks" function(req::HTTP.Request, id::Int)
        db   = current_db()
        rows = Tables.rowtable(DBInterface.execute(db, """
            SELECT a.id, a.exposure_id, a.q, a.intensity, a.prominence, a.sharpness,
                   'auto' AS source,
                   CASE WHEN c.q IS NOT NULL THEN 1 ELSE 0 END AS excluded
            FROM auto_peaks a
            LEFT JOIN peak_curations c
                ON c.exposure_id = a.exposure_id
               AND c.kind = 'exclude'
               AND ABS(c.q - a.q) <= MAX(1e-6, ABS(a.q) * 0.001)
            WHERE a.exposure_id = ?
            UNION ALL
            SELECT id, exposure_id, q, NULL AS intensity, NULL AS prominence, NULL AS sharpness,
                   'manual' AS source, 0 AS excluded
            FROM peak_curations
            WHERE exposure_id = ? AND kind = 'add'
            ORDER BY q
        """, [id, id]))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(rows_to_json(rows; bool_keys = (:excluded,))))
    end

    @post "/api/exposures/{id}/peaks" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)
        q    = Float64(body.q)

        apply_event!(db, req;
            kind        = "peak_added",
            entity_type = "exposure",
            entity_id   = id,
            payload     = Dict(:q => q))

        # Read back the just-created curation row.
        cur = first(Tables.rowtable(DBInterface.execute(db,
            """SELECT id FROM peak_curations
               WHERE exposure_id = ? AND kind = 'add' AND q = ?
               ORDER BY id DESC LIMIT 1""",
            [id, q])))
        peak_id = Int(cur.id)

        HTTP.Response(201, ["Content-Type" => "application/json"],
            JSON3.write(Dict(:id => peak_id, :exposure_id => id,
                             :q => q, :source => "manual", :excluded => false)))
    end

    @patch "/api/peaks/{id}" function(req::HTTP.Request, id::Int)
        db   = current_db()
        body = json(req)

        # Currently the only patchable field is `excluded`, and only on auto peaks.
        if !haskey(body, :excluded)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "missing field: excluded")))
        end

        # Check curation-add first: PATCH is not allowed on manual peaks
        # (auto peak ids and curation ids are in separate sequences and may
        # collide — always resolve 'add' curations before auto_peaks for PATCH).
        curation_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM peak_curations WHERE id = ? AND kind = 'add'", [id]))
        if !isempty(curation_rows)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "only auto peaks can be excluded; delete manual peaks instead")))
        end

        # Check auto_peaks
        auto_rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT a.id, a.exposure_id, a.q, a.intensity, a.prominence, a.sharpness,
                      'auto' AS source,
                      CASE WHEN c.q IS NOT NULL THEN 1 ELSE 0 END AS excluded
               FROM auto_peaks a
               LEFT JOIN peak_curations c
                   ON c.exposure_id = a.exposure_id
                  AND c.kind = 'exclude'
                  AND ABS(c.q - a.q) <= MAX(1e-6, ABS(a.q) * 0.001)
               WHERE a.id = ?""", [id]))
        if isempty(auto_rows)
            return HTTP.Response(404,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "peak not found")))
        end

        peak_row    = auto_rows[1]
        exposure_id = Int(peak_row.exposure_id)
        source      = "auto"

        excluded = Bool(body.excluded)
        peak_q   = Float64(peak_row.q)
        tol      = _peak_curation_tol(peak_q)

        if excluded
            # INSERT exclude curation — idempotent: only emit event if none exists
            existing = Tables.rowtable(DBInterface.execute(db,
                """SELECT id FROM peak_curations
                   WHERE exposure_id = ? AND kind = 'exclude'
                     AND ABS(q - ?) <= ?""",
                [exposure_id, peak_q, tol]))
            if isempty(existing)
                apply_event!(db, req;
                    kind        = "peak_excluded",
                    entity_type = "exposure",
                    entity_id   = exposure_id,
                    payload     = Dict(:q => peak_q, :auto_peak_id => id))
            end
        else
            # peak_unexcluded: find the prior peak_excluded event to record as undoes.
            prior = Tables.rowtable(DBInterface.execute(db, """
                SELECT id FROM user_actions
                WHERE action = 'peak_excluded'
                  AND entity_type = 'exposure' AND entity_id = ?
                  AND payload IS NOT NULL
                  AND ABS(json_extract(payload, '\$.q') - ?) <= MAX(1e-6, ABS(?) * 0.001)
                ORDER BY id DESC LIMIT 1
            """, [exposure_id, peak_q, peak_q]))
            undoes = isempty(prior) ? nothing : Int(prior[1].id)

            apply_event!(db, req;
                kind            = "peak_unexcluded",
                entity_type     = "exposure",
                entity_id       = exposure_id,
                payload         = Dict(:q => peak_q, :auto_peak_id => id),
                undoes_event_id = undoes)
        end

        # Return the updated row with fresh excluded state
        out = _effective_peak_row(db, id)
        out === nothing && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "peak not found after update")))

        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(row_to_json(out; bool_keys = (:excluded,))))
    end

    @delete "/api/peaks/{id}" function(req::HTTP.Request, id::Int)
        db = current_db()

        # Check peak_curations(kind='add') first — manual peaks are deletable.
        # (Auto peak ids and curation ids are independent sequences and can
        # collide; always prefer the curation interpretation for DELETE.)
        curation_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, exposure_id FROM peak_curations WHERE id = ? AND kind = 'add'",
            [id]))
        if !isempty(curation_rows)
            exposure_id = Int(curation_rows[1].exposure_id)

            # Delete the curation row and its associated index_peaks rows
            DBInterface.execute(db,
                "DELETE FROM index_peaks WHERE peak_id = ? AND peak_kind = 'curation'", [id])
            DBInterface.execute(db,
                "DELETE FROM peak_curations WHERE id = ?", [id])

            log_action!(db, req; action = "remove_peak",
                entity_type = "peak", entity_id = id)

            return HTTP.Response(204)
        end

        # If not a curation-add, check auto_peaks
        auto_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM auto_peaks WHERE id = ?", [id]))
        if !isempty(auto_rows)
            # Auto peaks cannot be hard-deleted; use PATCH { excluded: true } instead
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "auto peaks can only be excluded, not deleted; use PATCH { excluded: true } instead")))
        end

        HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "peak not found")))
    end

    @get "/api/peaks/{id}" function(req::HTTP.Request, id::Int)
        db  = current_db()
        row = _effective_peak_row(db, id)
        row === nothing && return HTTP.Response(404,
            ["Content-Type" => "application/json"],
            JSON3.write(Dict(:error => "peak not found")))
        HTTP.Response(200, ["Content-Type" => "application/json"],
            JSON3.write(row_to_json(row; bool_keys = (:excluded,))))
    end
end
