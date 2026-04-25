using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "sample_messages routes" begin
    tmp = mktempdir()
    db     = HimalayaUI.open_db(tmp)
    exp_id = HimalayaUI.create_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")

    with_test_server(db) do port, base
        # GET empty list
        r = HTTP.get("$base/api/samples/$s_id/messages")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test length(list) == 0

        # POST first message
        r = HTTP.post("$base/api/samples/$s_id/messages";
            body = JSON3.write(Dict(:body => "looks cubic to me")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 201
        msg = JSON3.read(String(r.body))
        @test msg.author == "alice"
        @test msg.body   == "looks cubic to me"
        @test msg.sample_id == s_id
        @test haskey(msg, :id)
        @test haskey(msg, :author_id)
        @test msg.author_id !== nothing
        @test haskey(msg, :created_at)

        # The FK should point at the user row created via get_or_create_user!
        user_rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM users WHERE username = ?", ["alice"]))
        @test !isempty(user_rows)
        @test msg.author_id == Int(user_rows[1].id)

        # POST a second message as a different user
        r = HTTP.post("$base/api/samples/$s_id/messages";
            body = JSON3.write(Dict(:body => "agreed — Im3m a=19.3 nm")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "bob"])
        @test r.status == 201

        # GET list — chronological (oldest first)
        r = HTTP.get("$base/api/samples/$s_id/messages")
        list = JSON3.read(String(r.body))
        @test length(list) == 2
        @test list[1].author == "alice"
        @test list[2].author == "bob"
        @test list[1].body   == "looks cubic to me"
        @test list[2].body   == "agreed — Im3m a=19.3 nm"

        # POST without X-Username → 401
        r = HTTP.post("$base/api/samples/$s_id/messages";
            body = JSON3.write(Dict(:body => "anon")),
            headers = ["Content-Type" => "application/json"],
            status_exception = false)
        @test r.status == 401

        # POST empty body → 400
        r = HTTP.post("$base/api/samples/$s_id/messages";
            body = JSON3.write(Dict(:body => "   ")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            status_exception = false)
        @test r.status == 400

        # GET for unknown sample → 200 empty (graceful)
        r = HTTP.get("$base/api/samples/99999/messages")
        @test r.status == 200
        @test length(JSON3.read(String(r.body))) == 0

        # POST to nonexistent sample → FK violation → non-2xx (FK enforcement is on)
        r = HTTP.post("$base/api/samples/99999/messages";
            body = JSON3.write(Dict(:body => "orphan message")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            status_exception = false)
        @test r.status >= 400

        # Deleting the user should null out `author` on subsequent reads (the LEFT JOIN
        # has nothing to resolve). `author_id` may or may not be preserved depending on
        # whether `PRAGMA foreign_keys = ON` is active — we only require that the UI's
        # display path (author name) sees null.
        DBInterface.execute(db, "DELETE FROM users WHERE username = ?", ["alice"])
        r = HTTP.get("$base/api/samples/$s_id/messages")
        list = JSON3.read(String(r.body))
        alice_msg = first(filter(m -> m.body == "looks cubic to me", list))
        @test alice_msg.author === nothing
    end
end
