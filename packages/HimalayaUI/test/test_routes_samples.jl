using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "samples routes" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/t", data_dir="/t/d",
                                             analysis_dir="/t/a")
    s1 = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    s2 = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D2", display_name="UX2")

    with_test_server(db) do port, base
        # List
        r = HTTP.get("$base/api/experiments/$exp_id/samples")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test length(list) == 2
        @test list[1].name == "D1"
        @test list[1].tags  == []

        # PATCH
        r = HTTP.patch("$base/api/samples/$s1";
            body = JSON3.write(Dict(:notes => "hello")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 200
        @test JSON3.read(String(r.body)).notes == "hello"

        # POST tag
        r = HTTP.post("$base/api/samples/$s1/tags";
            body = JSON3.write(Dict(:key => "lipid", :value => "DOPC")),
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 201
        tag = JSON3.read(String(r.body))
        @test tag.key   == "lipid"
        @test tag.value == "DOPC"
        tag_id = tag.id

        # Samples list now shows the tag
        r    = HTTP.get("$base/api/experiments/$exp_id/samples")
        list = JSON3.read(String(r.body))
        @test length(list[1].tags) == 1
        @test list[1].tags[1].key == "lipid"

        # DELETE tag
        r = HTTP.delete("$base/api/samples/$s1/tags/$tag_id";
            headers = ["X-Username" => "alice"])
        @test r.status == 204

        r    = HTTP.get("$base/api/experiments/$exp_id/samples")
        list = JSON3.read(String(r.body))
        @test list[1].tags == []
    end
end

@testset "PATCH /api/samples/:id rejects name (now immutable), accepts display_name" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/t3", data_dir="/t3/d",
                                             analysis_dir="/t3/a")
    sid = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX-immut")

    with_test_server(db) do port, base
        # :name is no longer in the allowlist — should return 400.
        r = HTTP.request("PATCH", "$base/api/samples/$sid",
            ["Content-Type" => "application/json",
             "X-Username"   => "alice"],
            JSON3.write(Dict(:name => "renamed")); status_exception=false)
        @test r.status == 400

        # :display_name is accepted; leading/trailing whitespace is trimmed.
        r2 = HTTP.request("PATCH", "$base/api/samples/$sid",
            ["Content-Type" => "application/json",
             "X-Username"   => "alice"],
            JSON3.write(Dict(:display_name => "  spaced  ")))
        @test r2.status == 200
        body = JSON3.read(String(r2.body))
        @test body[:display_name] == "spaced"
    end
end

@testset "POST /api/samples/:id/tags is idempotent under retry" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/t2", data_dir="/t2/d",
                                             analysis_dir="/t2/a")
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")

    with_test_server(db) do port, base
        op_id = "uuid-tag-retry-1"
        body  = JSON3.write(Dict(:key => "lipid", :value => "DOPC"))
        hdrs  = ["Content-Type"   => "application/json",
                 "X-Username"     => "alice",
                 "X-Client-Op-Id" => op_id]

        r1 = HTTP.post("$base/api/samples/$s_id/tags"; body=body, headers=hdrs)
        @test r1.status == 201
        body1 = String(copy(r1.body))

        r2 = HTTP.post("$base/api/samples/$s_id/tags"; body=body, headers=hdrs)
        @test r2.status == 201
        @test String(copy(r2.body)) == body1

        n_tags = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM sample_tags WHERE sample_id = ?", [s_id]))).c
        @test n_tags == 1

        n_events = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions WHERE action = 'add_tag' AND entity_type = 'sample'"))).c
        @test n_events == 1
    end
end
