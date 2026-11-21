using HTTP, JSON3, DBInterface, Tables, Oxygen

function register_picker_routes!()
    @get "/api/experiments/{eid}/sample-tags" function(req::HTTP.Request, eid::Int)
        db = current_db()
        # Per-(key, value) sample count joined through samples → experiment.
        # No need to validate the experiment id — an unknown id naturally
        # produces an empty list, which is the right answer for the picker.
        rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT t.key, t.value, COUNT(DISTINCT t.sample_id) AS count
               FROM sample_tags t
               JOIN samples s ON s.id = t.sample_id
               WHERE s.experiment_id = ?
               GROUP BY t.key, t.value
               ORDER BY t.key, t.value""", [eid]))

        out = [Dict(:key => String(r.key), :value => String(r.value), :count => Int(r.count)) for r in rows]
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
    end

    @get "/api/sample-tags" function(req::HTTP.Request)
        db = current_db()
        # Per-(key, value) sample count across the whole corpus, with no
        # experiment filter. COUNT(DISTINCT sample_id) gives the number of
        # samples carrying each (key, value) pair — used by the Manage-tags
        # modal to rank suggestions by frequency. proposeOrdering (scoping)
        # ranks by distinct-value count instead and ignores this field.
        rows = Tables.rowtable(DBInterface.execute(db,
            """SELECT key, value, COUNT(DISTINCT sample_id) AS count
               FROM sample_tags
               GROUP BY key, value
               ORDER BY key, value"""))

        out = [Dict(:key => String(r.key), :value => String(r.value), :count => Int(r.count)) for r in rows]
        HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(out))
    end
end
