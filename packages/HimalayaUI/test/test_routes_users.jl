using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "users routes" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)

    with_inproc_routes(db) do call
        # Empty list
        r = call("GET", "/api/users")
        @test r.status == 200
        @test JSON3.read(String(r.body)) == []

        # Create without first/last (legacy path)
        r = call("POST", "/api/users";
            headers = ["Content-Type" => "application/json"],
            body = Vector{UInt8}(JSON3.write(Dict(:username => "alice"))))
        @test r.status == 201
        created = JSON3.read(String(r.body))
        @test created.username == "alice"
        @test created.id == 1
        @test isnothing(created.first_name)
        @test isnothing(created.last_name)

        # Idempotent — second create returns existing
        r = call("POST", "/api/users";
            headers = ["Content-Type" => "application/json"],
            body = Vector{UInt8}(JSON3.write(Dict(:username => "alice"))))
        @test r.status == 200
        @test JSON3.read(String(r.body)).id == 1

        # Idempotent enrichment: a follow-up POST that supplies first/last fills in
        # the previously-NULL fields (but never overwrites non-null values).
        r = call("POST", "/api/users";
            headers = ["Content-Type" => "application/json"],
            body = Vector{UInt8}(JSON3.write(Dict(:username => "alice",
                                    :first_name => "Alice",
                                    :last_name  => "Smith"))))
        @test r.status == 200
        enriched = JSON3.read(String(r.body))
        @test enriched.id == 1
        @test enriched.first_name == "Alice"
        @test enriched.last_name  == "Smith"
        # Existing names are preserved on a subsequent POST with different names.
        r = call("POST", "/api/users";
            headers = ["Content-Type" => "application/json"],
            body = Vector{UInt8}(JSON3.write(Dict(:username => "alice",
                                    :first_name => "DIFFERENT",
                                    :last_name  => "OVERWRITE"))))
        @test r.status == 200
        preserved = JSON3.read(String(r.body))
        @test preserved.first_name == "Alice"
        @test preserved.last_name  == "Smith"

        # Create with first_name and last_name
        r = call("POST", "/api/users";
            headers = ["Content-Type" => "application/json"],
            body = Vector{UInt8}(JSON3.write(Dict(:username => "jwhc", :first_name => "Jonathan", :last_name => "Chen"))))
        @test r.status == 201
        created2 = JSON3.read(String(r.body))
        @test created2.username   == "jwhc"
        @test created2.first_name == "Jonathan"
        @test created2.last_name  == "Chen"

        # List returns first_name/last_name
        r = call("GET", "/api/users")
        users = JSON3.read(String(r.body))
        @test length(users) == 2
        jwhc = first(filter(u -> u.username == "jwhc", users))
        @test jwhc.first_name == "Jonathan"
        @test jwhc.last_name  == "Chen"

        # Reject invalid handles (spaces, @, hyphens, etc.)
        for bad in ("has space", "with-dash", "@withat", "with.dot")
            r = call("POST", "/api/users";
                headers = ["Content-Type" => "application/json"],
                body = Vector{UInt8}(JSON3.write(Dict(:username => bad))))
            @test r.status == 400
        end

        # Empty audit trail
        r = call("GET", "/api/users/alice/actions")
        @test r.status == 200
        @test JSON3.read(String(r.body)) == []

        # 404 for unknown user
        r = call("GET", "/api/users/nobody/actions")
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
    # note is now stored in the payload JSON column (via apply_event! delegation)
    @test JSON3.read(rows[1].payload).note == "hello"
end
