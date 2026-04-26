using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "PATCH /api/exposures/:id/status" begin
    db      = SQLite.DB()
    HimalayaUI.create_schema!(db)
    HimalayaUI.migrate_schema!(db)
    exp_id  = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
    samp_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
    eid     = HimalayaUI.create_exposure!(db; sample_id=samp_id)

    with_test_server(db) do port, base
        # accept
        r = HTTP.patch("$base/api/exposures/$eid/status";
            body    = JSON3.write(Dict(:status => "accepted")),
            headers = ["Content-Type" => "application/json"])
        @test r.status == 200
        row = Tables.rowtable(DBInterface.execute(db,
            "SELECT status FROM exposures WHERE id = ?", [eid]))[1]
        @test row.status == "accepted"

        # reject
        HTTP.patch("$base/api/exposures/$eid/status";
            body    = JSON3.write(Dict(:status => "rejected")),
            headers = ["Content-Type" => "application/json"])
        row2 = Tables.rowtable(DBInterface.execute(db,
            "SELECT status FROM exposures WHERE id = ?", [eid]))[1]
        @test row2.status == "rejected"

        # null (clear)
        HTTP.patch("$base/api/exposures/$eid/status";
            body    = JSON3.write(Dict(:status => nothing)),
            headers = ["Content-Type" => "application/json"])
        row3 = Tables.rowtable(DBInterface.execute(db,
            "SELECT status FROM exposures WHERE id = ?", [eid]))[1]
        @test row3.status === missing || row3.status === nothing

        # invalid value → 422
        r422 = HTTP.patch("$base/api/exposures/$eid/status";
            body    = JSON3.write(Dict(:status => "garbage")),
            headers = ["Content-Type" => "application/json"],
            status_exception = false)
        @test r422.status == 422

        # nonexistent exposure → 404
        r404 = HTTP.patch("$base/api/exposures/9999/status";
            body    = JSON3.write(Dict(:status => "accepted")),
            headers = ["Content-Type" => "application/json"],
            status_exception = false)
        @test r404.status == 404
    end
end
