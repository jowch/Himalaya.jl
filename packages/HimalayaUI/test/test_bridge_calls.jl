using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

if !isdefined(@__MODULE__, :with_test_server)
    include("test_http.jl")
end

function _bridge_call_fixture(dir)
    db = open_prepared_clone(dir)
    exp_id = HimalayaUI.create_experiment!(db; path = "/x", data_dir = "/x", analysis_dir = "/x")
    a = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "A")
    b = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "B")
    ea = HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = a)
    eb = HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = b)
    for (id, eid, phase, score) in ((10, ea, "Pn3m", 0.9), (11, ea, "Im3m", 0.5),
                                    (20, eb, "Pn3m", 0.8), (21, eb, "Im3m", 0.7))
        DBInterface.execute(db,
            "INSERT INTO indices (id, exposure_id, phase, basis, score) VALUES (?, ?, ?, 0.1, ?)",
            [id, eid, phase, score])
    end
    (; db, exp_id, a, b, ea, eb)
end

const _BC_HEADERS = ["Content-Type" => "application/json", "X-Username" => "alice"]
_bc_post(call, path, body) = call("POST", path; headers = _BC_HEADERS, body = Vector{UInt8}(JSON3.write(body)))
_bc_json(r) = JSON3.read(String(r.body))
_bc_tick(call, eid, iid) = _bc_post(call, "/api/exposures/$eid/assignment/members", Dict(:index_id => iid))
_bc_state(call, sid, state) = _bc_post(call, "/api/samples/$sid/call/state", Dict(:state => state))
_bc_bridge(call, a, b; keep = nothing) =
    _bc_post(call, "/api/samples/$a/bridge", keep === nothing ? Dict(:to_sample_id => b) : Dict(:to_sample_id => b, :keep_call => keep))

function _bc_listed(call, exp_id, sid)
    s = only(filter(s -> s.id == sid, collect(_bc_json(call("GET", "/api/samples?experiment_id=$exp_id")))))
    (state = String(s.assignment_state), phase = s.phase === nothing ? nothing : String(s.phase),
     phases = String.(collect(s.call_phases)))
end

_bc_members(db) = sort([(Int(r.exposure_id), Int(r.index_id)) for r in Tables.rowtable(
    DBInterface.execute(db, "SELECT exposure_id, index_id FROM assignment_members"))])

_bc_states(db) = Dict(Int(r.sample_id) => String(r.state) for r in Tables.rowtable(
    DBInterface.execute(db, "SELECT sample_id, state FROM sample_calls")))

@testset "bridged samples share one call" begin
    @testset "a phase ticked on one member's frame is the call of every member" begin
        mktempdir() do dir
            fx = _bridge_call_fixture(dir)
            with_inproc_routes(fx.db) do call
                @test _bc_bridge(call, fx.a, fx.b).status == 200
                @test _bc_tick(call, fx.ea, 10).status == 200
                for sid in (fx.a, fx.b)
                    @test _bc_listed(call, fx.exp_id, sid) == (state = "indexed", phase = "Pn3m", phases = ["Pn3m"])
                end
                body = _bc_json(call("GET", "/api/samples/$(fx.b)/call"))
                @test [(Int(m.exposure_id), Int(m.index_id)) for m in body.members] == [(fx.ea, 10)]
                @test String.(collect(body.call_phases)) == ["Pn3m"]
                @test _bc_json(call("GET", "/api/exposures/$(fx.eb)/assignment")).members == []
                exported = Dict(Int(s.sample_id) => s for s in _bc_json(call("GET", "/api/experiments/$(fx.exp_id)/export")))
                @test String(exported[fx.b].call_state) == "indexed"
                @test HimalayaUI.bridged_call_drift(fx.db) == Int[]
            end
        end
    end

    @testset "setting the call state on one member sets it on every member, through the log" begin
        mktempdir() do dir
            fx = _bridge_call_fixture(dir)
            with_inproc_routes(fx.db) do call
                _bc_bridge(call, fx.a, fx.b)
                _bc_tick(call, fx.ea, 10)
                _bc_tick(call, fx.eb, 21)
                @test _bc_state(call, fx.b, "form_factor").status == 200
                @test _bc_states(fx.db) == Dict(fx.a => "form_factor", fx.b => "form_factor")
                @test _bc_members(fx.db) == []
                @test _bc_listed(call, fx.exp_id, fx.a).state == "form_factor"
                logged = Tables.rowtable(DBInterface.execute(fx.db,
                    "SELECT entity_id FROM user_actions WHERE action = 'sample_call_set_state' ORDER BY id"))
                @test Set(Int(r.entity_id) for r in logged) == Set([fx.a, fx.b])
                @test _bc_tick(call, fx.eb, 20).status == 200
                @test _bc_states(fx.db) == Dict(fx.a => "indexed", fx.b => "indexed")
                @test HimalayaUI.bridged_call_drift(fx.db) == Int[]
            end
        end
    end

    @testset "bridging to an uncalled sample gives it the other side's call" begin
        mktempdir() do dir
            fx = _bridge_call_fixture(dir)
            with_inproc_routes(fx.db) do call
                _bc_state(call, fx.a, "form_factor")
                @test _bc_bridge(call, fx.b, fx.a).status == 200
                @test _bc_states(fx.db)[fx.b] == "form_factor"
                @test HimalayaUI.bridged_call_drift(fx.db) == Int[]
            end
        end
    end

    @testset "bridging two different calls is refused until the request says which to keep" begin
        mktempdir() do dir
            fx = _bridge_call_fixture(dir)
            with_inproc_routes(fx.db) do call
                _bc_tick(call, fx.ea, 10)
                _bc_tick(call, fx.eb, 21)
                r = _bc_bridge(call, fx.a, fx.b)
                @test r.status == 409
                calls = Dict(Int(c.sample_id) => String.(collect(c.call_phases)) for c in _bc_json(r).calls)
                @test calls == Dict(fx.a => ["Pn3m"], fx.b => ["Im3m"])
                @test isempty(_bc_json(call("GET", "/api/experiments/$(fx.exp_id)/bridges")))
                @test _bc_members(fx.db) == [(fx.ea, 10), (fx.eb, 21)]
                @test _bc_bridge(call, fx.a, fx.b; keep = 999).status == 400

                @test _bc_bridge(call, fx.a, fx.b; keep = fx.a).status == 200
                @test _bc_members(fx.db) == [(fx.ea, 10)]
                @test _bc_listed(call, fx.exp_id, fx.b).phases == ["Pn3m"]
                removed = Tables.rowtable(DBInterface.execute(fx.db,
                    "SELECT entity_id FROM user_actions WHERE action = 'assignment_remove'"))
                @test [Int(r.entity_id) for r in removed] == [fx.eb]
            end
        end
    end

    @testset "keeping both calls joins them, and both is refused when a side is not indexed" begin
        mktempdir() do dir
            fx = _bridge_call_fixture(dir)
            with_inproc_routes(fx.db) do call
                _bc_tick(call, fx.ea, 10)
                _bc_tick(call, fx.eb, 21)
                @test _bc_bridge(call, fx.a, fx.b; keep = "both").status == 200
                @test _bc_listed(call, fx.exp_id, fx.a).phases == ["Im3m", "Pn3m"]
                @test _bc_listed(call, fx.exp_id, fx.b).phases == ["Im3m", "Pn3m"]
            end
        end
        mktempdir() do dir
            fx = _bridge_call_fixture(dir)
            with_inproc_routes(fx.db) do call
                _bc_tick(call, fx.ea, 10)
                _bc_state(call, fx.b, "form_factor")
                @test _bc_bridge(call, fx.a, fx.b; keep = "both").status == 422
                @test _bc_bridge(call, fx.a, fx.b; keep = fx.b).status == 200
                @test _bc_states(fx.db) == Dict(fx.a => "form_factor", fx.b => "form_factor")
                @test _bc_members(fx.db) == []
            end
        end
    end

    @testset "leaving a bridge keeps the phases the sample's own frames support" begin
        mktempdir() do dir
            fx = _bridge_call_fixture(dir)
            with_inproc_routes(fx.db) do call
                _bc_bridge(call, fx.a, fx.b)
                _bc_tick(call, fx.ea, 10)
                _bc_tick(call, fx.eb, 21)
                @test call("DELETE", "/api/samples/$(fx.b)/bridge"; headers = _BC_HEADERS).status == 200
                @test _bc_listed(call, fx.exp_id, fx.a).phases == ["Pn3m"]
                @test _bc_listed(call, fx.exp_id, fx.b).phases == ["Im3m"]
            end
        end
    end

    @testset "the shared call survives a rebuild from the log" begin
        mktempdir() do dir
            fx = _bridge_call_fixture(dir)
            with_inproc_routes(fx.db) do call
                _bc_tick(call, fx.ea, 10)
                _bc_tick(call, fx.eb, 21)
                _bc_bridge(call, fx.a, fx.b; keep = fx.a)
                _bc_state(call, fx.a, "form_factor")
                _bc_state(call, fx.b, "indexed")
                _bc_tick(call, fx.eb, 20)
                before = (_bc_members(fx.db), _bc_states(fx.db))
                for sid in (fx.a, fx.b)
                    HimalayaUI.rebuild_views_from_log!(fx.db, sid; entity_type = "sample")
                end
                for eid in (fx.ea, fx.eb)
                    HimalayaUI.rebuild_views_from_log!(fx.db, eid; entity_type = "exposure")
                end
                @test (_bc_members(fx.db), _bc_states(fx.db)) == before
                @test HimalayaUI.bridged_call_drift(fx.db) == Int[]
            end
        end
    end

    @testset "bridged_call_drift names a group whose members' call states differ" begin
        mktempdir() do dir
            fx = _bridge_call_fixture(dir)
            with_inproc_routes(fx.db) do call
                _bc_bridge(call, fx.a, fx.b)
            end
            HimalayaUI._upsert_sample_call!(fx.db, fx.a, "form_factor")
            @test length(HimalayaUI.bridged_call_drift(fx.db)) == 1
        end
    end
end
