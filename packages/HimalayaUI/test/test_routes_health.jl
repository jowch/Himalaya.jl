@testset "GET /api/health" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)

    with_test_server(db) do port, base
        resp = HTTP.get("$base/api/health")
        @test resp.status == 200
        body = JSON3.read(String(resp.body))
        @test body.status == "ok"
    end
end

@testset "serve warns when it would run on one worker thread" begin
    @test_logs (:warn, r"--threads=4,1") @test HimalayaUI.warn_if_single_threaded(1)
    @test_logs @test !HimalayaUI.warn_if_single_threaded(4)
end
