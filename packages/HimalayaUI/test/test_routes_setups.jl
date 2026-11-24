using Test, HTTP, SQLite, DBInterface, Tables, HimalayaUI
import JSON
isdefined(@__MODULE__, :_write_setup!) || include("test_setups.jl")

@testset "PATCH /api/setups/:id" begin
    bt = mktempdir()
    d1, a1 = _write_setup!(joinpath(bt, "setup1"); stems = ["A_S0001_0_001"])
    d2, a2 = _write_setup!(joinpath(bt, "setup2"); stems = ["B_S0001_0_001"])
    db = fresh_db()
    eid = HimalayaUI.create_experiment!(db; name = "bt", path = bt, data_dir = d1, analysis_dir = a1)
    s1 = HimalayaUI.default_setup_id(db, eid)
    s2 = HimalayaUI.create_setup!(db; experiment_id = eid, data_dir = d2, analysis_dir = a2,
        position = 1, included = false)
    HimalayaUI.scan_and_group!(db, eid; analyze = false)
    sub = (pending = Channel{String}(64),)
    hdrs = ["Content-Type" => "application/json", "X-Username" => "alice"]
    patch(call, id, body) = call("PATCH", "/api/setups/$id"; headers = hdrs, body = Vector{UInt8}(body))
    with_inproc_routes(db) do call
        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end
        r = patch(call, s1, """{"name":"  Long pipe  "}""")
        @test r.status == 200
        @test JSON.parse(r.body).name == "Long pipe"
        @test HimalayaUI.setup_row(db, s1).display_name == "Long pipe"
        @test JSON.parse(patch(call, s1, """{"name":"setup1"}""").body).name == "setup1"
        @test ismissing(HimalayaUI.setup_row(db, s1).display_name)
        @test JSON.parse(patch(call, s1, """{"name":null}""").body).display_name === nothing

        @test patch(call, s1, """{"included":false}""").status == 400
        @test HimalayaUI.setup_row(db, s1).included == 1

        r = patch(call, s2, """{"included":true}""")
        @test r.status == 200
        @test JSON.parse(r.body).included === true
        t0 = time()
        while time() - t0 < 30
            n = only(Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ?", [s2]))).c
            n > 0 && !HimalayaUI.ingest_running(eid) && break
            sleep(0.05)
        end
        @test only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ?", [s2]))).c == 1
        @test patch(call, s1, """{"included":false}""").status == 200
        _write_setup!(joinpath(bt, "setup1"); stems = ["A_S0001_0_001", "A_S0002_0_001"])
        @test patch(call, s1, """{"included":true}""").status == 200
        t0 = time()
        while time() - t0 < 30
            n = only(Tables.rowtable(DBInterface.execute(db,
                "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ?", [s1]))).c
            n == 2 && !HimalayaUI.ingest_running(eid) && break
            sleep(0.05)
        end
        @test only(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS c FROM exposures WHERE setup_id = ?", [s1]))).c == 2

        for bad in ("""{"included":"no"}""", """{"name":3}""", """{"folder":"x"}""", "{}",
                    """{"name":"x","data_dir":"/y"}""")
            @test patch(call, s2, bad).status == 400
        end
        @test patch(call, 99999, """{"name":"x"}""").status == 404

        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT entity_type, entity_id, payload FROM user_actions WHERE action = 'setup_updated' ORDER BY id"))
        @test length(rows) == 6
        @test all(r -> r.entity_type == "setup", rows)
        @test JSON.parse(rows[1].payload).experiment_id == eid
        @test JSON.parse(rows[1].payload).name == "Long pipe"
        @test JSON.parse(rows[4].payload).included === true
    end
    frames = String[]
    while isready(sub.pending)
        push!(frames, take!(sub.pending))
    end
    @test count(f -> occursin("setup_updated", f), frames) == 6
end
