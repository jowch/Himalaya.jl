using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "samples routes" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/t", data_dir="/t/d",
                                             analysis_dir="/t/a")
    s1 = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    s2 = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D2")

    with_inproc_routes(db) do call
        # List
        r = call("GET", "/api/experiments/$exp_id/samples")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test length(list) == 2
        @test list[1].name == "D1"
        @test !haskey(list[1], :tags)

        # PATCH
        r = call("PATCH", "/api/samples/$s1";
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:notes => "hello"))))
        @test r.status == 200
        @test JSON3.read(String(r.body)).notes == "hello"

    end
end

@testset "GET /api/samples: per-sample phase-call rollup" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/t", data_dir="/t/d", analysis_dir="/t/a")

    # Sample A — INDEXED. Two exposures; the LOWER-id one is selected=1, so it is
    # the representative (proves selected wins over highest-id). Its assignment
    # carries Pn3m (score 0.9) + Lamellar (0.4) → dominant phase = Pn3m.
    sA  = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="A")
    eA1 = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=sA, filename="a1")
    HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=sA, filename="a2")   # higher id, NOT selected
    DBInterface.execute(db, "UPDATE exposures SET selected = 1 WHERE id = ?", [eA1])
    rPn = DBInterface.execute(db,
        "INSERT INTO indices (exposure_id, phase, basis, score) VALUES (?, 'Pn3m', 0.1, 0.9)", [eA1])
    idPn = Int(DBInterface.lastrowid(rPn))
    rLam = DBInterface.execute(db,
        "INSERT INTO indices (exposure_id, phase, basis, score) VALUES (?, 'Lamellar', 0.1, 0.4)", [eA1])
    idLam = Int(DBInterface.lastrowid(rLam))
    DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [eA1, idPn])
    DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [eA1, idLam])

    # Sample B — FORM FACTOR. The sample's call is declared form_factor.
    sB = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="B")
    eB = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=sB, filename="b1")
    DBInterface.execute(db, "INSERT INTO sample_calls (sample_id, state) VALUES (?, 'form_factor')", [sB])

    # Sample C — UNINDEXED. Representative exposure, no assignment row.
    sC = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="C")
    eC = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=sC, filename="c1")

    with_inproc_routes(db) do call
        list   = JSON3.read(String(call("GET", "/api/samples").body))
        byname = Dict(String(s.name) => s for s in list)

        @test byname["A"].assignment_state == "indexed"
        @test byname["A"].phase == "Pn3m"        # top-scored member of the sample's call
        @test collect(byname["A"].call_phases) == ["Lamellar", "Pn3m"]
        @test byname["B"].assignment_state == "form_factor"
        @test byname["B"].phase === nothing      # form factor carries no phase
        @test byname["C"].assignment_state == "indexed"  # default (no members)
        @test byname["C"].phase === nothing      # unindexed
        @test byname["A"].representative_exposure_id == eA1
        @test byname["B"].representative_exposure_id == eB
        @test byname["C"].representative_exposure_id == eC
    end
end

@testset "PATCH /api/samples/:id accepts name, leading/trailing whitespace trimmed" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/t3", data_dir="/t3/d",
                                             analysis_dir="/t3/a")
    sid = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")

    with_inproc_routes(db) do call
        # :name is accepted; leading/trailing whitespace is trimmed.
        r2 = call("PATCH", "/api/samples/$sid";
            headers = ["Content-Type" => "application/json",
                       "X-Username"   => "alice"],
            body = Vector{UInt8}(JSON3.write(Dict(:name => "  spaced  "))))
        @test r2.status == 200
        body = JSON3.read(String(r2.body))
        @test body[:name] == "spaced"
    end
end

@testset "corpus samples route" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)

    e1 = HimalayaUI.init_experiment!(db; name="E1", path="/e1",
        data_dir="/e1/d", analysis_dir="/e1/a")
    e2 = HimalayaUI.init_experiment!(db; name="E2", path="/e2",
        data_dir="/e2/d", analysis_dir="/e2/a")

    DBInterface.execute(db, "UPDATE setups SET q_units = 'nm-1' WHERE experiment_id = ?", [e1])

    s1 = HimalayaUI.create_sample!(db; experiment_id=e1, name="A1")
    s2 = HimalayaUI.create_sample!(db; experiment_id=e1, name="A2")
    s3 = HimalayaUI.create_sample!(db; experiment_id=e2, name="B1")

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
        @test by_name["B1"].q_units == "A-1"

        @test !any(haskey(s, :tags) for s in all)

        # ?experiment_id= filter narrows to one experiment.
        r = call("GET", "/api/samples?experiment_id=$e1")
        @test r.status == 200
        filtered = JSON3.read(String(r.body))
        @test length(filtered) == 2
        @test Set(s.name for s in filtered) == Set(["A1", "A2"])
        @test JSON3.write.(filtered) == JSON3.write.([s for s in all if s.experiment_id == e1])

        # Nonexistent experiment id → empty array (SQL gives this for free).
        r = call("GET", "/api/samples?experiment_id=999999")
        @test r.status == 200
        @test JSON3.read(String(r.body)) == []

        # Malformed experiment_id → 400, not a silent ignore.
        r = call("GET", "/api/samples?experiment_id=abc")
        @test r.status == 400
    end
end

@testset "corpus samples route reads q_units from the sample's setup" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)

    e_typed = HimalayaUI.create_experiment!(db; name="ET", path="/et",
        data_dir="/et/d", analysis_dir="/et/a", q_units="nm^-1")
    e_blob = HimalayaUI.init_experiment!(db; name="EB", path="/eb",
        data_dir="/eb/d", analysis_dir="/eb/a")
    e_both = HimalayaUI.create_experiment!(db; name="EX", path="/ex",
        data_dir="/ex/d", analysis_dir="/ex/a", q_units="A^-1")

    HimalayaUI.create_sample!(db; experiment_id=e_typed, name="T1")
    HimalayaUI.create_sample!(db; experiment_id=e_blob,  name="L1")
    HimalayaUI.create_sample!(db; experiment_id=e_both,  name="X1")

    with_inproc_routes(db) do call
        r = call("GET", "/api/samples")
        @test r.status == 200
        by_name = Dict(String(s.name) => s for s in JSON3.read(String(r.body)))
        @test by_name["T1"].q_units == "nm^-1"   
        @test by_name["L1"].q_units == "A-1"
        @test by_name["X1"].q_units == "A^-1"     
    end
end

@testset "corpus samples route tolerates NULL experiment_id" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)

    # `samples.experiment_id` is a nullable FK. A sample with no owning
    # experiment (a SQL NULL) must not 500 the corpus route — its q_units
    # falls back to the default rather than throwing in `Int(...)`.
    DBInterface.execute(db,
        "INSERT INTO samples (experiment_id, name)
         VALUES (NULL, 'orphan')")

    with_inproc_routes(db) do call
        r = call("GET", "/api/samples")
        @test r.status == 200
        rows = JSON3.read(String(r.body))
        @test length(rows) == 1
        @test rows[1].q_units == "A-1"
    end
end

@testset "representative frame: shared cross-view cases" begin
    cases = JSON3.read(read(joinpath(@__DIR__, "..", "frontend", "test", "fixtures",
                                     "representative-exposure.json"), String))
    for c in cases
        @testset "$(c.case)" begin
            @test HimalayaUI.representative_exposure_id(c.frames) == c.representative
        end
    end
end

@testset "GET /api/samples: the representative is never a dropped frame" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/t", data_dir="/t/d", analysis_dir="/t/a")
    frames(name, n) = (s = HimalayaUI.create_sample!(db; experiment_id=exp_id, name=name);
        [HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s, filename="$(name)$(i)") for i in 1:n])
    drop(e) = DBInterface.execute(db, "UPDATE exposures SET status = 'rejected' WHERE id = ?", [e])

    newest = frames("newest_dropped", 3); drop(newest[3])
    picked = frames("selected_dropped", 3); drop(picked[1])
    DBInterface.execute(db, "UPDATE exposures SET selected = 1 WHERE id = ?", [picked[1]])
    gone = frames("all_dropped", 2); drop.(gone)

    with_inproc_routes(db) do call
        byname = Dict(String(s.name) => s for s in JSON3.read(String(call("GET", "/api/samples").body)))
        @test byname["newest_dropped"].representative_exposure_id == newest[2]
        @test byname["selected_dropped"].representative_exposure_id == picked[3]
        @test byname["all_dropped"].representative_exposure_id === nothing
        @test byname["all_dropped"].assignment_state == "indexed"
        @test byname["all_dropped"].phase === nothing
    end
end
