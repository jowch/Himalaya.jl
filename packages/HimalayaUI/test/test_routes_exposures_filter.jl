using Test, HTTP, JSON3, SQLite, DBInterface

@testset "GET /api/samples/:id/exposures exclude_rejected" begin
    db      = SQLite.DB()
    HimalayaUI.create_schema!(db)
    HimalayaUI.migrate_schema!(db)
    exp_id  = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
    samp_id = HimalayaUI.create_sample!(db; experiment_id=exp_id)
    HimalayaUI.create_exposure!(db; sample_id=samp_id, filename="good.dat")
    HimalayaUI.create_exposure!(db; sample_id=samp_id, filename="bad.dat", status="rejected")

    with_test_server(db) do port, base
        # without filter — both returned
        r_all = HTTP.get("$base/api/samples/$samp_id/exposures")
        all_exps = JSON3.read(String(r_all.body))
        @test length(all_exps) == 2

        # with filter — rejected excluded
        r_fil = HTTP.get("$base/api/samples/$samp_id/exposures?exclude_rejected=true")
        fil_exps = JSON3.read(String(r_fil.body))
        @test length(fil_exps) == 1
        @test String(fil_exps[1].filename) == "good.dat"
    end
end
