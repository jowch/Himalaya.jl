using Test, HTTP, JSON3, SQLite, DBInterface, Tables

@testset "exposures routes" begin
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
    e1 = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="example_tot",
        exposure_time=15.0, transmitted_counts=7365.0, incident_monitor=1.5007, counting_time=15.8)
    e2 = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename="other")
    prp = joinpath(tmp, "with_prp.prp")
    write(prp, "Time this file was written: 26 Apr 2026 18:20:27\nBeam energy=9000.5 eV\nPipe length=1700 mm\nDetector=Pilatus 1M\n")
    s2_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D2")
    e3 = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s2_id, filename="with_prp",
        prp_path=prp, timestamp="2026-04-26T18:20:27")
    e4 = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s2_id, filename="prp_gone",
        prp_path=joinpath(tmp, "missing.prp"))

    with_inproc_routes(db) do call
        # List
        r = call("GET", "/api/samples/$s_id/exposures")
        @test r.status == 200
        list = JSON3.read(String(r.body))
        @test length(list) == 2
        @test list[1].filename == "example_tot"
        @test list[1].selected === true
        @test !haskey(list[1], :tags)
        @test list[1].sources  == []
        @test list[1].transmitted_counts ≈ 7365.0
        @test list[1].incident_monitor ≈ 1.5007
        @test list[1].counting_time ≈ 15.8
        @test list[1].exposure_time ≈ 15.0
        @test list[2].transmitted_counts === nothing
        @test list[1].beam_energy_ev === nothing
        @test list[1].pipe_length_m === nothing
        @test list[1].detector === nothing

        r = call("GET", "/api/samples/$s2_id/exposures")
        @test r.status == 200
        withprp = JSON3.read(String(r.body))
        @test withprp[1].id == e3
        @test withprp[1].timestamp == "2026-04-26T18:20:27"
        @test withprp[1].beam_energy_ev ≈ 9000.5
        @test withprp[1].pipe_length_m ≈ 1.7
        @test withprp[1].detector == "Pilatus 1M"
        @test withprp[2].id == e4
        @test withprp[2].beam_energy_ev === nothing
        @test withprp[2].pipe_length_m === nothing
        @test withprp[2].detector === nothing

        write(prp, "Beam energy=1.0 eV\nPipe length=1 mm\nDetector=rewritten\n")
        write(joinpath(tmp, "missing.prp"), "Beam energy=11000 eV\nDetector=late\n")
        again = JSON3.read(String(call("GET", "/api/samples/$s2_id/exposures").body))
        @test again[1].beam_energy_ev ≈ 9000.5
        @test again[1].detector == "Pilatus 1M"
        @test again[2].beam_energy_ev ≈ 11000
        @test again[2].detector == "late"

        # Analyze single exposure
        r = call("POST", "/api/exposures/$e1/analyze";
            headers = ["X-Username" => "alice"])
        @test r.status == 200
        @test JSON3.read(String(r.body)).analyzed === true

        # 404 for nonexistent
        r = call("POST", "/api/exposures/9999/analyze";
            headers = ["X-Username" => "alice"])
        @test r.status == 404
    end
end

@testset "prp_details: a failed parse is not cached" begin
    path = joinpath(mktempdir(), "flaky.prp")
    write(path, "Beam energy=8000 eV\n")
    @test (@test_logs (:warn, r"could not read") match_mode=:any HimalayaUI.prp_details(path; parse = _ -> error("unreadable")))[:beam_energy_ev] === nothing
    @test HimalayaUI.prp_details(path)[:beam_energy_ev] ≈ 8000
    @test HimalayaUI.prp_details(path; parse = _ -> error("not called"))[:beam_energy_ev] ≈ 8000
end
