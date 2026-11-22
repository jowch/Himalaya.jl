using Test, Dates, SQLite, DBInterface, Tables, JSON3, HimalayaUI
isdefined(@__MODULE__, :_write_setup!) || include("test_setups.jl")

_bt_rows(db, sql, args...) = Tables.rowtable(DBInterface.execute(db, sql, collect(args)))

function _wait_ingest(db, eid; timeout = 60)
    t0 = time()
    while time() - t0 < timeout
        st = only(_bt_rows(db, "SELECT ingest_status FROM experiments WHERE id = ?", eid)).ingest_status
        isequal(st, "scanning") || HimalayaUI.ingest_running(eid) || return st
        sleep(0.05)
    end
    error("ingest did not finish")
end

_post(call, body) = call("POST", "/api/experiments";
    headers = ["Content-Type" => "application/json"], body = Vector{UInt8}(JSON3.write(body)))

@testset "beamtime workspaces" begin
    @testset "a rescan adds a new setup folder to the beamtime" begin
        bt = mktempdir()
        d1, a1 = _write_setup!(joinpath(bt, "setup1"); stems = ["A_S0001_0_001"])
        db = fresh_db()
        eid = HimalayaUI.create_experiment!(db; name = "bt", path = realpath(bt),
            data_dir = realpath(d1), analysis_dir = realpath(a1))
        @test isempty(HimalayaUI.new_setup_roots(db, eid))
        HimalayaUI.scan_and_group!(db, eid; analyze = false)

        _write_setup!(joinpath(bt, "setup2"); stems = ["B_S0001_0_001"])
        @test HimalayaUI.new_setup_roots(db, eid) == [realpath(joinpath(bt, "setup2"))]
        @test HimalayaUI.ingest_has_work(db, eid)

        added = HimalayaUI.add_new_setups!(db, eid)
        @test length(added) == 1
        st = HimalayaUI.setup_row(db, only(added))
        @test st.folder == "setup2"
        @test st.position == 1
        @test st.included == 1
        @test st.data_dir == realpath(joinpath(bt, "setup2", "data"))
        @test isempty(HimalayaUI.add_new_setups!(db, eid))

        HimalayaUI.scan_and_group!(db, eid; analyze = false)
        @test _count(db, "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ?", st.id) == 1
    end

    @testset "an unticked setup is remembered and not re-added" begin
        bt = mktempdir()
        d1, a1 = _write_setup!(joinpath(bt, "setup1"); stems = ["A_S0001_0_001"])
        d2, a2 = _write_setup!(joinpath(bt, "setup2"); stems = ["B_S0001_0_001"])
        db = fresh_db()
        eid = HimalayaUI.create_experiment!(db; name = "bt", path = realpath(bt),
            data_dir = realpath(d1), analysis_dir = realpath(a1))
        HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = realpath(d2),
            analysis_dir = realpath(a2), position = 1, included = false)
        @test isempty(HimalayaUI.new_setup_roots(db, eid))
        HimalayaUI.scan_and_group!(db, eid; analyze = false)
        @test _count(db, "SELECT COUNT(*) AS c FROM exposures WHERE experiment_id = ?", eid) == 1
    end

    @testset "a one-setup experiment never adopts subfolders" begin
        bt = mktempdir()
        root = joinpath(bt, "setup1")
        d1, a1 = _write_setup!(root; stems = ["A_S0001_0_001"])
        _write_setup!(joinpath(d1, "nested"); stems = ["N_S0001_0_001"])
        db = fresh_db()
        eid = HimalayaUI.create_experiment!(db; name = "s1", path = realpath(root),
            data_dir = realpath(d1), analysis_dir = realpath(a1))
        @test isempty(HimalayaUI.new_setup_roots(db, eid))
    end

    @testset "a folder another experiment uses is skipped" begin
        bt = mktempdir()
        d1, a1 = _write_setup!(joinpath(bt, "setup1"); stems = ["A_S0001_0_001"])
        d2, a2 = _write_setup!(joinpath(bt, "setup2"); stems = ["B_S0001_0_001"])
        db = fresh_db()
        HimalayaUI.create_experiment!(db; name = "old", path = realpath(joinpath(bt, "setup2")),
            data_dir = realpath(d2), analysis_dir = realpath(a2))
        eid = HimalayaUI.create_experiment!(db; name = "bt", path = realpath(bt),
            data_dir = realpath(d1), analysis_dir = realpath(a1))
        @test isempty(HimalayaUI.new_setup_roots(db, eid))
    end

    @testset "an unreachable beamtime folder adds nothing" begin
        db = fresh_db()
        eid = HimalayaUI.create_experiment!(db; name = "gone", path = "/nonexistent/beamtime")
        @test isempty(HimalayaUI.new_setup_roots(db, eid))
    end

    @testset "POST /api/experiments with setups makes one beamtime" begin
        bt = mktempdir()
        d1, a1 = _write_setup!(joinpath(bt, "setup1"); stems = ["A_S0001_0_001"])
        d2, a2 = _write_setup!(joinpath(bt, "setup2"); stems = ["B_S0001_0_001"])
        db = fresh_db()
        with_inproc_routes(db) do call
            r = _post(call, Dict(:name => "  Nov beamtime ", :path => bt, :setups => [
                Dict(:root => joinpath(bt, "setup1"), :data_dir => d1, :analysis_dir => a1,
                     :patterns => Dict(:image => "{name}.tif"),
                     :geometry => Dict(:flight_path_m => 1.7, :flight_path_m_source => "setup")),
                Dict(:root => joinpath(bt, "setup2"), :data_dir => d2, :analysis_dir => a2,
                     :included => false)]))
            @test r.status == 202
            eid = Int(JSON3.read(String(r.body)).id)
            @test _wait_ingest(db, eid) == "complete"
            x = only(_bt_rows(db, "SELECT name, path FROM experiments WHERE id = ?", eid))
            @test x.name == "Nov beamtime"
            @test x.path == realpath(bt)
            sts = HimalayaUI.setups_of(db, eid)
            @test [s.folder for s in sts] == ["setup1", "setup2"]
            @test [s.included for s in sts] == [1, 0]
            @test sts[1].image_pattern == "{name}.tif"
            @test sts[1].flight_path_m_source == "setup"
            @test _count(db, "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ?", sts[2].id) == 0
            @test _count(db, "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ?", sts[1].id) == 1

            body = JSON3.read(String(call("GET", "/api/experiments/$eid").body))
            @test [s.included for s in body.setups] == [true, false]

            dup = _post(call, Dict(:path => bt, :setups => [Dict(:data_dir => d2)]))
            @test dup.status == 409
            @test _count(db, "SELECT COUNT(*) AS c FROM experiments") == 1

            none = _post(call, Dict(:path => mktempdir(), :setups => [Dict(:data_dir => mktempdir(), :included => false)]))
            @test none.status == 400
            @test _post(call, Dict(:path => bt, :setups => [])).status == 400
            @test _count(db, "SELECT COUNT(*) AS c FROM experiments") == 1
        end
    end

    @testset "merge folds per-setup experiments into one beamtime" begin
        bt = mktempdir()
        d1, a1 = _write_setup!(joinpath(bt, "setup1"); stems = ["A_S0001_0_001", "A_S0002_0_001"])
        d2, a2 = _write_setup!(joinpath(bt, "setup2"); stems = ["B_S0001_0_001"])
        db = fresh_db()
        e2 = HimalayaUI.create_experiment!(db; name = "two", path = realpath(joinpath(bt, "setup2")),
            data_dir = realpath(d2), analysis_dir = realpath(a2))
        e1 = HimalayaUI.create_experiment!(db; name = "one", path = realpath(joinpath(bt, "setup1")),
            data_dir = realpath(d1), analysis_dir = realpath(a1))
        HimalayaUI.scan_and_group!(db, e1; analyze = false)
        HimalayaUI.scan_and_group!(db, e2; analyze = false)
        samples_before = sort([r.id for r in _bt_rows(db, "SELECT id FROM samples")])
        exposures_before = sort([r.id for r in _bt_rows(db, "SELECT id FROM exposures")])

        io = IOBuffer()
        plan = HimalayaUI.cli_merge(["two", "one"]; io, db)
        out = String(take!(io))
        @test occursin("Dry run", out)
        @test plan.into == e2
        @test plan.path == realpath(bt)
        @test plan.name == basename(realpath(bt))
        @test [s.folder for s in plan.setups] == ["setup1", "setup2"]
        @test _count(db, "SELECT COUNT(*) AS c FROM experiments") == 2

        HimalayaUI.cli_merge([string(e2), string(e1), "--name", "Nov", "--apply"]; io = IOBuffer(), db)
        @test _count(db, "SELECT COUNT(*) AS c FROM experiments") == 1
        x = only(_bt_rows(db, "SELECT id, name, path FROM experiments"))
        @test (x.id, x.name, x.path) == (e2, "Nov", realpath(bt))
        @test [s.folder for s in HimalayaUI.setups_of(db, e2)] == ["setup1", "setup2"]
        @test sort([r.id for r in _bt_rows(db, "SELECT id FROM samples WHERE experiment_id = ?", e2)]) == samples_before
        @test sort([r.id for r in _bt_rows(db, "SELECT id FROM exposures WHERE experiment_id = ?", e2)]) == exposures_before
        @test _count(db, "SELECT COUNT(*) AS c FROM loads WHERE experiment_id != ?", e2) == 0
        @test isempty(HimalayaUI.new_setup_roots(db, e2))
        again = HimalayaUI.scan_and_group!(db, e2; analyze = false)
        @test again.added_exposures == 0

        @test_throws ErrorException HimalayaUI.merge_plan(db, [e2])
        @test_throws ErrorException HimalayaUI.merge_plan(db, [e2, e2])
    end

    @testset "merge refuses two setups with the same folder name" begin
        db = fresh_db()
        a = HimalayaUI.create_experiment!(db; name = "a", path = "/x/day1/setup1",
            data_dir = "/x/day1/setup1/data", analysis_dir = "/x/day1/setup1/analysis")
        b = HimalayaUI.create_experiment!(db; name = "b", path = "/x/day2/setup1",
            data_dir = "/x/day2/setup1/data", analysis_dir = "/x/day2/setup1/analysis")
        @test_throws ErrorException HimalayaUI.merge_plan(db, [a, b])
    end
end
