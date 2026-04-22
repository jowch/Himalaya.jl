using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "users routes" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)

    with_test_server(db) do port, base
        # Empty list
        r = HTTP.get("$base/api/users")
        @test r.status == 200
        @test JSON3.read(String(r.body)) == []

        # Create
        r = HTTP.post("$base/api/users";
            body = JSON3.write(Dict(:username => "alice")),
            headers = ["Content-Type" => "application/json"])
        @test r.status == 201
        created = JSON3.read(String(r.body))
        @test created.username == "alice"
        @test created.id == 1

        # Idempotent — second create returns existing
        r = HTTP.post("$base/api/users";
            body = JSON3.write(Dict(:username => "alice")),
            headers = ["Content-Type" => "application/json"])
        @test r.status == 200
        @test JSON3.read(String(r.body)).id == 1

        # List now has alice
        r = HTTP.get("$base/api/users")
        users = JSON3.read(String(r.body))
        @test length(users) == 1
        @test users[1].username == "alice"

        # Empty audit trail
        r = HTTP.get("$base/api/users/alice/actions")
        @test r.status == 200
        @test JSON3.read(String(r.body)) == []

        # 404 for unknown user
        r = HTTP.get("$base/api/users/nobody/actions"; status_exception = false)
        @test r.status == 404
    end
end

@testset "log_action!" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    HimalayaUI.get_or_create_user!(db, "alice")

    req = HTTP.Request("POST", "/dummy", ["X-Username" => "alice"])
    HimalayaUI.log_action!(db, req;
        action      = "test",
        entity_type = "sample",
        entity_id   = 42,
        note        = "hello")

    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT * FROM user_actions ORDER BY id DESC LIMIT 1"))
    @test length(rows) == 1
    @test rows[1].action      == "test"
    @test rows[1].entity_type == "sample"
    @test rows[1].entity_id   == 42
    @test rows[1].note        == "hello"
end
