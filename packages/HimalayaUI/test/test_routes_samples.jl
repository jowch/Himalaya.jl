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

@testset "POST /api/samples/tags/batch inserts N tags in one transaction" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tb", data_dir="/tb/d",
                                             analysis_dir="/tb/a")
    sA = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="DA")
    sB = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="DB")
    sC = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="DC")

    with_test_server(db) do port, base
        body = JSON3.write(Dict(
            :key    => "temperature",
            :source => "scoping",
            :tags   => [Dict(:sample_id => sA, :value => "25C"),
                        Dict(:sample_id => sB, :value => "30C"),
                        Dict(:sample_id => sC, :value => "35C")]))
        r = HTTP.post("$base/api/samples/tags/batch"; body = body,
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 201

        created = JSON3.read(String(r.body))
        @test length(created) == 3
        # Response order matches request order; each entry is a full SampleTag
        # plus the explicit sample_id (no URL param to make it implicit).
        @test created[1].sample_id == sA
        @test created[1].key       == "temperature"
        @test created[1].value     == "25C"
        @test created[1].source    == "scoping"
        @test created[3].sample_id == sC
        @test created[3].value     == "35C"
        @test all(t -> t.id > 0, created)

        # All three tags persisted.
        n_tags = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM sample_tags WHERE key = 'temperature'"))).c
        @test n_tags == 3

        # N add_tag events written — the route reuses the existing kind.
        n_events = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions
             WHERE action = 'add_tag' AND entity_type = 'sample'"))).c
        @test n_events == 3

        # Samples list reflects the new tags.
        list = JSON3.read(String(HTTP.get("$base/api/experiments/$exp_id/samples").body))
        @test length(list[1].tags) == 1
        @test list[1].tags[1].source == "scoping"
    end
end

@testset "POST /api/samples/tags/batch source defaults to manual" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tbd", data_dir="/tbd/d",
                                             analysis_dir="/tbd/a")
    sA = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="DA")

    with_test_server(db) do port, base
        body = JSON3.write(Dict(
            :key  => "lipid",
            :tags => [Dict(:sample_id => sA, :value => "DOPC")]))
        r = HTTP.post("$base/api/samples/tags/batch"; body = body,
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"])
        @test r.status == 201
        @test JSON3.read(String(r.body))[1].source == "manual"
    end
end

@testset "POST /api/samples/tags/batch rolls the whole batch back on a mid-batch failure" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    # Mirror production (open_db): FK enforcement makes a nonexistent
    # sample_id fail the INSERT — the mid-batch failure under test.
    DBInterface.execute(db, "PRAGMA foreign_keys = ON")
    exp_id = HimalayaUI.init_experiment!(db; path="/tbr", data_dir="/tbr/d",
                                             analysis_dir="/tbr/a")
    sA = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="DA")
    sB = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="DB")
    bad_sample_id = sB + 9999  # no such sample — FK violation mid-loop

    with_test_server(db) do port, base
        body = JSON3.write(Dict(
            :key  => "temperature",
            :tags => [Dict(:sample_id => sA,            :value => "25C"),
                      Dict(:sample_id => bad_sample_id, :value => "30C"),
                      Dict(:sample_id => sB,            :value => "35C")]))
        r = HTTP.request("POST", "$base/api/samples/tags/batch",
            ["Content-Type" => "application/json", "X-Username" => "alice"],
            body; status_exception = false)
        @test r.status >= 400

        # Atomic boundary: the valid sA insert before the bad entry must NOT
        # survive — no half-confirmed recipe.
        n_tags = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM sample_tags WHERE key = 'temperature'"))).c
        @test n_tags == 0

        n_events = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions WHERE action = 'add_tag'"))).c
        @test n_events == 0
    end
end

@testset "POST /api/samples/tags/batch is idempotent under retry" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tbi", data_dir="/tbi/d",
                                             analysis_dir="/tbi/a")
    sA = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="DA")
    sB = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="DB")

    with_test_server(db) do port, base
        op_id = "uuid-batch-tags-retry-1"
        body  = JSON3.write(Dict(
            :key    => "temperature",
            :source => "scoping",
            :tags   => [Dict(:sample_id => sA, :value => "25C"),
                        Dict(:sample_id => sB, :value => "30C")]))
        hdrs  = ["Content-Type"   => "application/json",
                 "X-Username"     => "alice",
                 "X-Client-Op-Id" => op_id]

        r1 = HTTP.post("$base/api/samples/tags/batch"; body = body, headers = hdrs)
        @test r1.status == 201
        body1 = String(copy(r1.body))

        r2 = HTTP.post("$base/api/samples/tags/batch"; body = body, headers = hdrs)
        @test r2.status == 201
        @test String(copy(r2.body)) == body1

        # The retry replayed the cached response — no duplicate rows.
        n_tags = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM sample_tags WHERE key = 'temperature'"))).c
        @test n_tags == 2

        n_events = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions WHERE action = 'add_tag'"))).c
        @test n_events == 2
    end
end

@testset "POST /api/samples/tags/batch rejects malformed requests" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tbv", data_dir="/tbv/d",
                                             analysis_dir="/tbv/a")
    sA = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="DA")

    with_test_server(db) do port, base
        post(b) = HTTP.request("POST", "$base/api/samples/tags/batch",
            ["Content-Type" => "application/json", "X-Username" => "alice"],
            JSON3.write(b); status_exception = false)

        # Missing key.
        @test post(Dict(:tags => [Dict(:sample_id => sA, :value => "x")])).status == 400
        # Missing tags.
        @test post(Dict(:key => "temperature")).status == 400
        # Empty tags array.
        @test post(Dict(:key => "temperature", :tags => [])).status == 400
        # Tag entry missing value.
        @test post(Dict(:key => "temperature",
                        :tags => [Dict(:sample_id => sA)])).status == 400
        # Tag entry missing sample_id.
        @test post(Dict(:key => "temperature",
                        :tags => [Dict(:value => "x")])).status == 400

        # Wrong-typed fields must be a clean 400, never an unguarded-
        # conversion 500 (test_route_validation_routing.jl's invariant).
        # Non-array tags.
        @test post(Dict(:key => "temperature", :tags => 5)).status == 400
        # Non-string key.
        @test post(Dict(:key => 5,
                        :tags => [Dict(:sample_id => sA, :value => "x")])).status == 400
        # Non-integer sample_id (string and fractional).
        @test post(Dict(:key => "temperature",
                        :tags => [Dict(:sample_id => "abc", :value => "x")])).status == 400
        @test post(Dict(:key => "temperature",
                        :tags => [Dict(:sample_id => 5.5, :value => "x")])).status == 400
        # Non-string value.
        @test post(Dict(:key => "temperature",
                        :tags => [Dict(:sample_id => sA, :value => 5)])).status == 400
        # Non-string source.
        @test post(Dict(:key => "temperature", :source => 5,
                        :tags => [Dict(:sample_id => sA, :value => "x")])).status == 400

        # Duplicate sample_id in one batch is rejected: every event shares one
        # client_op_id, so a repeated sample_id would collide on
        # idx_events_unique_op and leave a tag row with no durable add_tag
        # event. Reject it up front instead.
        @test post(Dict(:key => "temperature",
                        :tags => [Dict(:sample_id => sA, :value => "25C"),
                                  Dict(:sample_id => sA, :value => "30C")])).status == 400

        # None of the rejected requests wrote anything.
        n_tags = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM sample_tags"))).c
        @test n_tags == 0
    end
end
