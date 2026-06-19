using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "samples routes" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/t", data_dir="/t/d",
                                             analysis_dir="/t/a")
    s1 = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")
    s2 = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D2", display_name="UX2")

    with_inproc_routes(db) do call
        # List
        r = call("GET", "/api/experiments/$exp_id/samples")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test length(list) == 2
        @test list[1].name == "D1"
        @test list[1].tags  == []

        # PATCH
        r = call("PATCH", "/api/samples/$s1";
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:notes => "hello"))))
        @test r.status == 200
        @test JSON3.read(String(r.body)).notes == "hello"

        # POST tag
        r = call("POST", "/api/samples/$s1/tags";
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:key => "lipid", :value => "DOPC"))))
        @test r.status == 201
        tag = JSON3.read(String(r.body))
        @test tag.key   == "lipid"
        @test tag.value == "DOPC"
        tag_id = tag.id

        # Samples list now shows the tag
        r    = call("GET", "/api/experiments/$exp_id/samples")
        list = JSON3.read(String(r.body))
        @test length(list[1].tags) == 1
        @test list[1].tags[1].key == "lipid"

        # DELETE tag
        r = call("DELETE", "/api/samples/$s1/tags/$tag_id";
            headers = ["X-Username" => "alice"])
        @test r.status == 204

        r    = call("GET", "/api/experiments/$exp_id/samples")
        list = JSON3.read(String(r.body))
        @test list[1].tags == []
    end
end

@testset "GET /api/samples: per-sample indexing rollup (representative exposure)" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/t", data_dir="/t/d", analysis_dir="/t/a")

    # Sample A — INDEXED. Two exposures; the LOWER-id one is selected=1, so it is
    # the representative (proves selected wins over highest-id). Its assignment
    # carries Pn3m (score 0.9) + Lamellar (0.4) → dominant phase = Pn3m.
    sA  = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="A")
    eA1 = HimalayaUI.create_exposure!(db; sample_id=sA, filename="a1")
    HimalayaUI.create_exposure!(db; sample_id=sA, filename="a2")   # higher id, NOT selected
    DBInterface.execute(db, "UPDATE exposures SET selected = 1 WHERE id = ?", [eA1])
    rPn = DBInterface.execute(db,
        "INSERT INTO indices (exposure_id, phase, basis, score) VALUES (?, 'Pn3m', 0.1, 0.9)", [eA1])
    idPn = Int(DBInterface.lastrowid(rPn))
    rLam = DBInterface.execute(db,
        "INSERT INTO indices (exposure_id, phase, basis, score) VALUES (?, 'Lamellar', 0.1, 0.4)", [eA1])
    idLam = Int(DBInterface.lastrowid(rLam))
    DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [eA1, idPn])
    DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [eA1, idLam])

    # Sample B — FORM FACTOR. Representative exposure declared form_factor.
    sB = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="B")
    eB = HimalayaUI.create_exposure!(db; sample_id=sB, filename="b1")
    DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'form_factor')", [eB])

    # Sample C — UNINDEXED. Representative exposure, no assignment row.
    sC = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="C")
    HimalayaUI.create_exposure!(db; sample_id=sC, filename="c1")

    with_inproc_routes(db) do call
        list   = JSON3.read(String(call("GET", "/api/samples").body))
        byname = Dict(String(s.name) => s for s in list)

        @test byname["A"].assignment_state == "indexed"
        @test byname["A"].phase == "Pn3m"        # top-scored member of the SELECTED rep
        @test byname["B"].assignment_state == "form_factor"
        @test byname["B"].phase === nothing      # form factor carries no phase
        @test byname["C"].assignment_state == "indexed"  # default (no members)
        @test byname["C"].phase === nothing      # unindexed
    end
end

@testset "PATCH /api/samples/:id rejects name (now immutable), accepts display_name" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/t3", data_dir="/t3/d",
                                             analysis_dir="/t3/a")
    sid = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX-immut")

    with_inproc_routes(db) do call
        # :name is no longer in the allowlist — should return 400.
        r = call("PATCH", "/api/samples/$sid";
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:name => "renamed"))))
        @test r.status == 400

        # :display_name is accepted; leading/trailing whitespace is trimmed.
        r2 = call("PATCH", "/api/samples/$sid";
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:display_name => "  spaced  "))))
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

    with_inproc_routes(db) do call
        op_id = "uuid-tag-retry-1"
        body  = JSON3.write(Dict(:key => "lipid", :value => "DOPC"))
        hdrs  = ["Content-Type"   => "application/json",
                 "X-Username"     => "alice",
                 "X-Client-Op-Id" => op_id]

        r1 = call("POST", "/api/samples/$s_id/tags"; headers=hdrs, body=Vector{UInt8}(body))
        @test r1.status == 201
        body1 = String(copy(r1.body))

        r2 = call("POST", "/api/samples/$s_id/tags"; headers=hdrs, body=Vector{UInt8}(body))
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

@testset "POST /api/samples/:id/tags rejects a duplicate key with 409" begin
    # Single-valued-key invariant (TAG-DEDUP-MODEL): a second add of an
    # already-present key is a clean 409, not an unhandled 500 from the unique
    # index — and writes nothing. A distinct op_id (not a retry) so the
    # idempotency cache doesn't replay the first response.
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tdup409", data_dir="/tdup409/d",
                                             analysis_dir="/tdup409/a")
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")

    with_inproc_routes(db) do call
        hdrs = ["Content-Type" => "application/json", "X-Username" => "alice"]
        r1 = call("POST", "/api/samples/$s_id/tags";
            headers=hdrs,
            body=Vector{UInt8}(JSON3.write(Dict(:key => "dose", :value => "10"))))
        @test r1.status == 201

        # Same key, different value, fresh op_id → collision.
        r2 = call("POST", "/api/samples/$s_id/tags";
            headers=hdrs,
            body=Vector{UInt8}(JSON3.write(Dict(:key => "dose", :value => "20"))))
        @test r2.status == 409

        # Still exactly one dose row, original value untouched; no add_tag twin.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT value FROM sample_tags WHERE sample_id = ? AND key = 'dose'", [s_id]))
        @test length(rows) == 1
        @test String(rows[1].value) == "10"
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

    with_inproc_routes(db) do call
        body = JSON3.write(Dict(
            :key    => "temperature",
            :source => "scoping",
            :tags   => [Dict(:sample_id => sA, :value => "25C"),
                        Dict(:sample_id => sB, :value => "30C"),
                        Dict(:sample_id => sC, :value => "35C")]))
        r = call("POST", "/api/samples/tags/batch";
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            body = Vector{UInt8}(body))
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
        list = JSON3.read(String(call("GET", "/api/experiments/$exp_id/samples").body))
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

    with_inproc_routes(db) do call
        body = JSON3.write(Dict(
            :key  => "lipid",
            :tags => [Dict(:sample_id => sA, :value => "DOPC")]))
        r = call("POST", "/api/samples/tags/batch";
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            body = Vector{UInt8}(body))
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

    with_inproc_routes(db) do call
        body = JSON3.write(Dict(
            :key  => "temperature",
            :tags => [Dict(:sample_id => sA,            :value => "25C"),
                      Dict(:sample_id => bad_sample_id, :value => "30C"),
                      Dict(:sample_id => sB,            :value => "35C")]))
        r = call("POST", "/api/samples/tags/batch";
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
            body = Vector{UInt8}(body))
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

    with_inproc_routes(db) do call
        op_id = "uuid-batch-tags-retry-1"
        body  = JSON3.write(Dict(
            :key    => "temperature",
            :source => "scoping",
            :tags   => [Dict(:sample_id => sA, :value => "25C"),
                        Dict(:sample_id => sB, :value => "30C")]))
        hdrs  = ["Content-Type"   => "application/json",
                 "X-Username"     => "alice",
                 "X-Client-Op-Id" => op_id]

        r1 = call("POST", "/api/samples/tags/batch"; headers=hdrs, body=Vector{UInt8}(body))
        @test r1.status == 201
        body1 = String(copy(r1.body))

        r2 = call("POST", "/api/samples/tags/batch"; headers=hdrs, body=Vector{UInt8}(body))
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

    with_inproc_routes(db) do call
        hdrs = ["Content-Type" => "application/json", "X-Username" => "alice"]
        post(b) = call("POST", "/api/samples/tags/batch";
            headers = hdrs,
            body = Vector{UInt8}(JSON3.write(b)))

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

@testset "PATCH /api/samples/:id/tags/:tag_id edits in place" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tpe", data_dir="/tpe/d",
                                             analysis_dir="/tpe/a")
    sid = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S1")
    # Seed a tag directly.
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, ?)",
        [sid, "dose", "10", "manual"])
    T = Int(DBInterface.lastrowid(DBInterface.execute(db, "SELECT last_insert_rowid()")))

    with_inproc_routes(db) do call
        r = call("PATCH", "/api/samples/$sid/tags/$T";
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:value => "12"))))
        @test r.status == 200
        body = JSON3.read(String(r.body))
        @test body.id == T && body.value == "12" && body.key == "dose"
        # Row updated in place, id stable.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, value FROM sample_tags WHERE id = ?", [T]))
        @test length(rows) == 1 && rows[1].value == "12"
    end
end

@testset "PATCH /api/samples/:id/tags/:tag_id rejects key-collision with 409" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tp409", data_dir="/tp409/d",
                                             analysis_dir="/tp409/a")
    sid = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S1")
    # Seed two tags: dose=10 (A) and temp=25 (B).
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, ?)",
        [sid, "dose", "10", "manual"])
    A = Int(DBInterface.lastrowid(DBInterface.execute(db, "SELECT last_insert_rowid()")))
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, ?)",
        [sid, "temp", "25", "manual"])
    B = Int(DBInterface.lastrowid(DBInterface.execute(db, "SELECT last_insert_rowid()")))

    with_inproc_routes(db) do call
        # Editing B's key to "dose" collides with A.
        r = call("PATCH", "/api/samples/$sid/tags/$B";
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:key => "dose"))))
        @test r.status == 409
    end
end

@testset "PATCH value edit succeeds on a sample with a same-key duplicate" begin
    # The single-valued-key collision check must guard only key *changes*. A
    # value-only edit (or a no-op key write, which the client sends) on a sample
    # that ALREADY carries two same-key tags — the duplicate this modal exists
    # to resolve — must NOT 409: it doesn't change the key set.
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tpdup", data_dir="/tpdup/d",
                                             analysis_dir="/tpdup/a")
    sid = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S1")
    # Two byte-identical-key tags: dose=10 (manual) and dose=10 (scoping).
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, ?)",
        [sid, "dose", "10", "manual"])
    A = Int(DBInterface.lastrowid(DBInterface.execute(db, "SELECT last_insert_rowid()")))
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, ?)",
        [sid, "dose", "10", "scoping"])
    B = Int(DBInterface.lastrowid(DBInterface.execute(db, "SELECT last_insert_rowid()")))

    with_inproc_routes(db) do call
        # Value-only edit of A (key omitted): allowed despite B sharing the key.
        r = call("PATCH", "/api/samples/$sid/tags/$A";
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:value => "5"))))
        @test r.status == 200
        # No-op key write alongside a value edit (what the modal actually sends).
        r2 = call("PATCH", "/api/samples/$sid/tags/$B";
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:key => "dose", :value => "20"))))
        @test r2.status == 200
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT value FROM sample_tags WHERE id = ? ", [A]))
        @test String(rows[1].value) == "5"
    end
end

@testset "PATCH /api/samples/:id/tags/:tag_id returns 404 for non-matching tag" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tp404", data_dir="/tp404/d",
                                             analysis_dir="/tp404/a")
    sid = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S1")

    with_inproc_routes(db) do call
        r = call("PATCH", "/api/samples/$sid/tags/999999";
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:value => "x"))))
        @test r.status == 404
    end
end

@testset "edit_tag is a non-view log event (user_actions row, no view write)" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tplog", data_dir="/tplog/d",
                                             analysis_dir="/tplog/a")
    sid = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="S1")
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, ?)",
        [sid, "dose", "10", "manual"])
    T = Int(DBInterface.lastrowid(DBInterface.execute(db, "SELECT last_insert_rowid()")))

    with_inproc_routes(db) do call
        n0 = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions"))).c
        # Capture the list of view tables before the PATCH.
        view_counts0 = map(["indices","auto_peaks","peak_curations","assignment_members"]) do tbl
            only(Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM $tbl"))).c
        end

        call("PATCH", "/api/samples/$sid/tags/$T";
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:value => "13"))))

        n1 = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions"))).c
        @test n1 == n0 + 1

        # Confirm the event kind is correct.
        row = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT action FROM user_actions ORDER BY id DESC LIMIT 1")))
        @test row.action == "edit_tag"

        # No view tables were written.
        view_counts1 = map(["indices","auto_peaks","peak_curations","assignment_members"]) do tbl
            only(Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM $tbl"))).c
        end
        @test view_counts1 == view_counts0
    end
end

# ---------------------------------------------------------------------------
# Batch upsert tests (Slice 5 / Task 10)
# ---------------------------------------------------------------------------

@testset "POST /api/samples/tags/batch upserts: existing key updates in place, provenance preserved" begin
    # Sample S already has dose=10 (id X, source manual). A scoping batch
    # with key=dose value="12" must update the row in place — no twin — and
    # must NOT change the source (manual provenance preserved).
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tbu1", data_dir="/tbu1/d",
                                             analysis_dir="/tbu1/a")
    S = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="SU1")
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, ?)",
        [S, "dose", "10", "manual"])
    X = Int(DBInterface.lastrowid(DBInterface.execute(db, "SELECT last_insert_rowid()")))

    with_inproc_routes(db) do call
        n0 = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions"))).c

        r = call("POST", "/api/samples/tags/batch";
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:key => "dose", :source => "scoping",
                                   :tags => [Dict(:sample_id => S, :value => "12")]))))
        @test r.status == 201

        # Exactly one dose row remains, same id, updated value.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, value, source FROM sample_tags
             WHERE sample_id = ? AND key = 'dose'", [S]))
        @test length(rows) == 1
        @test rows[1].id    == X
        @test String(rows[1].value) == "12"
        # Provenance untouched — scoping must not overwrite manual.
        @test String(rows[1].source) == "manual"

        # An edit_tag event was written (not add_tag).
        n1 = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions"))).c
        @test n1 == n0 + 1
        last_action = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT action FROM user_actions ORDER BY id DESC LIMIT 1")))
        @test last_action.action == "edit_tag"

        # Response carries one entry with the resolved id, new value, actual source.
        body = JSON3.read(String(r.body))
        @test length(body) == 1
        @test body[1].id     == X
        @test body[1].value  == "12"
        @test body[1].source == "manual"
    end
end

@testset "POST /api/samples/tags/batch upserts: unchanged value is a no-op" begin
    # Sample S has dose=12 (from a prior batch/insert). Re-scoping with the
    # same value must write nothing and emit no event.
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tbu2", data_dir="/tbu2/d",
                                             analysis_dir="/tbu2/a")
    S = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="SU2")
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, ?)",
        [S, "dose", "12", "scoping"])
    X = Int(DBInterface.lastrowid(DBInterface.execute(db, "SELECT last_insert_rowid()")))

    with_inproc_routes(db) do call
        n0 = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions"))).c

        r = call("POST", "/api/samples/tags/batch";
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:key => "dose", :source => "scoping",
                                   :tags => [Dict(:sample_id => S, :value => "12")]))))
        @test r.status == 201

        # No new event: value was unchanged.
        n1 = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions"))).c
        @test n1 == n0

        # Still exactly one row, same id, same value.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, value FROM sample_tags
             WHERE sample_id = ? AND key = 'dose'", [S]))
        @test length(rows) == 1
        @test rows[1].id == X
        @test String(rows[1].value) == "12"

        # Response still carries one entry per input row.
        body = JSON3.read(String(r.body))
        @test length(body) == 1
        @test body[1].id == X
    end
end

@testset "POST /api/samples/tags/batch upserts: new key inserts" begin
    # A key the sample does not yet have: should INSERT, emit add_tag.
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tbu3", data_dir="/tbu3/d",
                                             analysis_dir="/tbu3/a")
    S = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="SU3")

    with_inproc_routes(db) do call
        n0 = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions"))).c

        r = call("POST", "/api/samples/tags/batch";
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:key => "lipid", :source => "scoping",
                                   :tags => [Dict(:sample_id => S, :value => "DOPC")]))))
        @test r.status == 201

        # New row inserted, source matches the batch.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, value, source FROM sample_tags
             WHERE sample_id = ? AND key = 'lipid'", [S]))
        @test length(rows) == 1
        @test String(rows[1].value)  == "DOPC"
        @test String(rows[1].source) == "scoping"

        # add_tag event written.
        n1 = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions"))).c
        @test n1 == n0 + 1
        last_action = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT action FROM user_actions ORDER BY id DESC LIMIT 1")))
        @test last_action.action == "add_tag"

        body = JSON3.read(String(r.body))
        @test length(body) == 1
        @test body[1].id > 0
        @test body[1].source == "scoping"
    end
end

@testset "POST /api/samples/tags/batch upsert: idempotency replay writes nothing new" begin
    # Same client_op_id twice: the second call returns the cached response and
    # must produce no new rows or events beyond the first call.
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tbu4", data_dir="/tbu4/d",
                                             analysis_dir="/tbu4/a")
    S = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="SU4")
    # Seed an existing tag to exercise the update-branch idempotency.
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, ?)",
        [S, "dose", "5", "manual"])

    with_inproc_routes(db) do call
        op_id = "uuid-batch-upsert-idem-1"
        body_str = JSON3.write(Dict(:key => "dose", :source => "scoping",
                                    :tags => [Dict(:sample_id => S, :value => "99")]))
        hdrs = ["Content-Type"   => "application/json",
                "X-Username"     => "alice",
                "X-Client-Op-Id" => op_id]

        r1 = call("POST", "/api/samples/tags/batch"; headers=hdrs, body=Vector{UInt8}(body_str))
        @test r1.status == 201
        body1 = String(copy(r1.body))

        n_tags0 = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM sample_tags WHERE sample_id = ? AND key = 'dose'",
            [S]))).c
        n_events0 = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions"))).c

        r2 = call("POST", "/api/samples/tags/batch"; headers=hdrs, body=Vector{UInt8}(body_str))
        @test r2.status == 201
        # Cached response replayed verbatim.
        @test String(copy(r2.body)) == body1

        # No new rows, no new events on the replay.
        n_tags1 = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM sample_tags WHERE sample_id = ? AND key = 'dose'",
            [S]))).c
        n_events1 = only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM user_actions"))).c
        @test n_tags1   == n_tags0
        @test n_events1 == n_events0
    end
end

@testset "corpus samples route" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)

    e1 = HimalayaUI.init_experiment!(db; name="E1", path="/e1",
        data_dir="/e1/d", analysis_dir="/e1/a")
    e2 = HimalayaUI.init_experiment!(db; name="E2", path="/e2",
        data_dir="/e2/d", analysis_dir="/e2/a")

    # Distinct q_units per experiment, written straight into the config blob.
    # e1 gets an explicit nm-1; e2 is left config-less so it falls back to A-1.
    DBInterface.execute(db,
        "UPDATE experiments SET config = ? WHERE id = ?",
        ["[beamline]\nq_units = \"nm-1\"\n", e1])

    s1 = HimalayaUI.create_sample!(db; experiment_id=e1, name="A1", display_name="UX-A1")
    s2 = HimalayaUI.create_sample!(db; experiment_id=e1, name="A2", display_name="UX-A2")
    s3 = HimalayaUI.create_sample!(db; experiment_id=e2, name="B1", display_name="UX-B1")

    # One tag on s1, so the projection's bundled `tags` array is exercised.
    DBInterface.execute(db,
        "INSERT INTO sample_tags (sample_id, key, value, source)
         VALUES (?, 'lipid', 'DOPC', 'manual')", [s1])

    with_inproc_routes(db) do call
        # Full corpus: every sample across both experiments.
        r = call("GET", "/api/samples")
        @test r.status == 200
        all = JSON3.read(String(r.body))
        @test length(all) == 3
        @test Set(s.name for s in all) == Set(["A1", "A2", "B1"])

        by_name = Dict(String(s.name) => s for s in all)

        # q_units sourced from each sample's owning experiment.
        @test by_name["A1"].q_units == "nm-1"
        @test by_name["A2"].q_units == "nm-1"
        @test by_name["B1"].q_units == "A-1"   # e2 has no config → default

        # tags bundled in the projection.
        @test length(by_name["A1"].tags) == 1
        @test by_name["A1"].tags[1].key   == "lipid"
        @test by_name["A1"].tags[1].value == "DOPC"
        @test by_name["A2"].tags == []

        # ?experiment_id= filter narrows to one experiment.
        r = call("GET", "/api/samples?experiment_id=$e1")
        @test r.status == 200
        filtered = JSON3.read(String(r.body))
        @test length(filtered) == 2
        @test Set(s.name for s in filtered) == Set(["A1", "A2"])

        # Nonexistent experiment id → empty array (SQL gives this for free).
        r = call("GET", "/api/samples?experiment_id=999999")
        @test r.status == 200
        @test JSON3.read(String(r.body)) == []

        # Malformed experiment_id → 400, not a silent ignore.
        r = call("GET", "/api/samples?experiment_id=abc")
        @test r.status == 400
    end
end

@testset "corpus samples route tolerates NULL experiment_id" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)

    # `samples.experiment_id` is a nullable FK. A sample with no owning
    # experiment (a SQL NULL) must not 500 the corpus route — its q_units
    # falls back to the default rather than throwing in `Int(...)`.
    DBInterface.execute(db,
        "INSERT INTO samples (experiment_id, name, display_name)
         VALUES (NULL, 'orphan', 'Orphan')")

    with_inproc_routes(db) do call
        r = call("GET", "/api/samples")
        @test r.status == 200
        rows = JSON3.read(String(r.body))
        @test length(rows) == 1
        @test rows[1].q_units == "A-1"
        @test rows[1].tags == []
    end
end

@testset "POST /api/samples/:id/tags — source defaults to manual, accepts an explicit value, rejects non-string" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tsrc", data_dir="/tsrc/d",
                                             analysis_dir="/tsrc/a")
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")

    with_inproc_routes(db) do call
        hdrs = ["Content-Type" => "application/json", "X-Username" => "alice"]

        # No source field → defaults to 'manual'.
        r = call("POST", "/api/samples/$s_id/tags";
            headers = hdrs,
            body = Vector{UInt8}(JSON3.write(Dict(:key => "lipid", :value => "DOPC"))))
        @test r.status == 201
        @test JSON3.read(String(r.body)).source == "manual"

        # Explicit source is honored end to end.
        r = call("POST", "/api/samples/$s_id/tags";
            headers = hdrs,
            body = Vector{UInt8}(JSON3.write(Dict(:key => "temp", :value => "25C",
                                    :source => "scoping"))))
        @test r.status == 201
        @test JSON3.read(String(r.body)).source == "scoping"

        # Non-string source → 400, nothing written.
        r = call("POST", "/api/samples/$s_id/tags";
            headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:key => "x", :value => "y", :source => 42))))
        @test r.status == 400

        # Exactly the two successful inserts persisted, with the right source.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT key, source FROM sample_tags WHERE sample_id = ? ORDER BY id",
            [s_id]))
        @test length(rows) == 2
        @test rows[1].source == "manual"
        @test rows[2].source == "scoping"
    end
end
