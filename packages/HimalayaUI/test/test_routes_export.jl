using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "export routes" begin
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id,
        name="D1")
    e_id   = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
    for (label, state) in (("Buffer", "confirmed"), ("maybe", "suggested"), ("37C", "confirmed"), ("gone", "removed"))
        DBInterface.execute(db, "INSERT INTO sample_labels (sample_id, label, state) VALUES (?, ?, ?)", [s_id, label, state])
    end

    with_inproc_routes(db) do call
        # JSON
        r = call("GET", "/api/experiments/$exp_id/export?format=json")
        @test r.status == 200
        @test occursin("application/json", HTTP.header(r, "Content-Type"))
        body = JSON3.read(String(r.body))
        @test length(body) == 1
        s = body[1]
        @test s.name == "D1"
        @test length(s.exposures) == 1
        @test s.exposures[1].filename == "example_tot"
        @test length(s.exposures[1].indices) >= 1
        @test collect(s.labels) == ["Buffer", "37C"]
        @test s.call_state == "uncalled"
        @test s.phase === nothing
        @test !s.exposures[1].called
        @test !haskey(s, :tags)

        # CSV
        r = call("GET", "/api/experiments/$exp_id/export?format=csv")
        @test r.status == 200
        @test occursin("text/csv", HTTP.header(r, "Content-Type"))
        csv_body = String(r.body)
        @test startswith(csv_body, "sample_name,exposure_filename,call_state,call_phases,phases")
        @test occursin("D1,example_tot,uncalled,,", csv_body)
        @test strip(last(split(first(split(csv_body, "\n")), ","))) == "labels"
        @test endswith(strip(split(csv_body, "\n")[2]), ",Buffer;37C")

        # Default = json
        r = call("GET", "/api/experiments/$exp_id/export")
        @test occursin("application/json", HTTP.header(r, "Content-Type"))

        # Invalid format
        r = call("GET", "/api/experiments/$exp_id/export?format=xml")
        @test r.status == 400
    end

    auto_ids = [Int(r.id) for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM indices WHERE exposure_id = ?", [e_id]))]
    called_id = Int(DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO indices (exposure_id, phase, basis, score, lattice_d, kind) VALUES (?, 'Im3m', 0.05, 0.1, 177.7, 'speculative')",
        [e_id])))
    DBInterface.execute(db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [e_id, called_id])
    DBInterface.execute(db, "INSERT INTO sample_calls (sample_id, state) VALUES (?, 'indexed')", [s_id])
    @test !isempty(Tables.rowtable(DBInterface.execute(db,
        "SELECT 1 FROM index_groups WHERE exposure_id = ? AND active = 1", [e_id])))

    with_inproc_routes(db) do call
        s = JSON3.read(String(call("GET", "/api/experiments/$exp_id/export").body))[1]
        @test s.call_state == "indexed"
        @test s.phase == "Im3m"
        @test collect(s.call_phases) == ["Im3m"]
        ex = s.exposures[1]
        @test ex.called
        @test [ix.phase for ix in ex.indices] == ["Im3m"]
        @test ex.indices[1].lattice_d == 177.7
        @test !haskey(ex, :active_group_kind)

        lines = split(strip(String(call("GET", "/api/experiments/$exp_id/export?format=csv").body)), "\n")
        @test startswith(lines[1], "sample_name,exposure_filename,call_state,call_phases,phases,")
        @test startswith(lines[2], "D1,example_tot,indexed,Im3m,Im3m,177.7,")
    end

    DBInterface.execute(db, "UPDATE sample_calls SET state = 'form_factor' WHERE sample_id = ?", [s_id])
    DBInterface.execute(db, "DELETE FROM assignment_members WHERE exposure_id = ?", [e_id])

    with_inproc_routes(db) do call
        s = JSON3.read(String(call("GET", "/api/experiments/$exp_id/export").body))[1]
        @test s.call_state == "form_factor"
        @test s.phase === nothing
        @test isempty(s.call_phases)
        ex = s.exposures[1]
        @test !ex.called
        @test sort([ix.phase for ix in ex.indices]) ==
              sort([String(r.phase) for r in Tables.rowtable(DBInterface.execute(db,
                  "SELECT phase FROM indices WHERE exposure_id = ?", [e_id]))])
        @test length(ex.indices) == length(auto_ids) + 1
    end
end
