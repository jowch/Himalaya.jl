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
