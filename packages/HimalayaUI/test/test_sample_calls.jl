using Test
using HimalayaUI
using SQLite, DBInterface, Tables
using JSON3, HTTP

function _call_fixture(dir; frames = 2)
    db = open_prepared_clone(dir)
    exp_id = HimalayaUI.create_experiment!(db; path = "/x", data_dir = "/x", analysis_dir = "/x")
    s_id = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = "S")
    eids = [HimalayaUI.create_exposure!(db; experiment_id = exp_id, sample_id = s_id) for _ in 1:frames]
    (; db, exp_id, s_id, eids)
end

_index!(db, id, eid, phase, score) = DBInterface.execute(db,
    "INSERT INTO indices (id, exposure_id, phase, basis, score) VALUES (?, ?, ?, 0.1, ?)",
    [id, eid, phase, score])

_json(r) = JSON3.read(String(r.body))
_post(call, path, body) = call("POST", path;
    headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
    body = Vector{UInt8}(JSON3.write(body)))
_listed(call, exp_id, s_id) = only(filter(s -> s.id == s_id,
    collect(_json(call("GET", "/api/samples?experiment_id=$exp_id")))))
_members(db) = sort([(Int(r.exposure_id), Int(r.index_id)) for r in Tables.rowtable(
    DBInterface.execute(db, "SELECT exposure_id, index_id FROM assignment_members"))])

if isdefined(@__MODULE__, :with_test_server)
    @testset "a phase call survives dropping the frame it was made on (#350)" begin
        mktempdir() do dir
            fx = _call_fixture(dir)
            e1, e2 = fx.eids
            _index!(fx.db, 10, e1, "Pn3m", 0.9)
            with_inproc_routes(fx.db) do call
                @test _post(call, "/api/exposures/$e1/assignment/members", Dict(:index_id => 10)).status == 200
                @test _listed(call, fx.exp_id, fx.s_id).phase == "Pn3m"
                r = call("PATCH", "/api/exposures/$e1/status";
                    headers = ["Content-Type" => "application/json", "X-Username" => "alice"],
                    body = Vector{UInt8}(JSON3.write(Dict(:status => "rejected"))))
                @test r.status == 200
                s = _listed(call, fx.exp_id, fx.s_id)
                @test s.assignment_state == "indexed"
                @test s.phase == "Pn3m"
                @test _json(call("GET", "/api/exposures/$e2/assignment")).state == "indexed"
            end
        end
    end

    @testset "a sample's call is the union of the phases picked on each frame" begin
        mktempdir() do dir
            fx = _call_fixture(dir)
            e1, e2 = fx.eids
            _index!(fx.db, 10, e1, "Pn3m", 0.5)
            _index!(fx.db, 20, e2, "HII", 0.9)
            with_inproc_routes(fx.db) do call
                _post(call, "/api/exposures/$e1/assignment/members", Dict(:index_id => 10))
                _post(call, "/api/exposures/$e2/assignment/members", Dict(:index_id => 20))
                @test _members(fx.db) == [(e1, 10), (e2, 20)]
                s = _listed(call, fx.exp_id, fx.s_id)
                @test s.phase == "HII"
                @test collect(s.call_phases) == ["HII", "Pn3m"]
                @test collect(_json(call("GET", "/api/exposures/$e1/assignment")).members) == [10]
            end
        end
    end

    @testset "setting the sample's call state applies to every frame" begin
        mktempdir() do dir
            fx = _call_fixture(dir)
            e1, e2 = fx.eids
            _index!(fx.db, 10, e1, "Pn3m", 0.5)
            _index!(fx.db, 20, e2, "HII", 0.9)
            with_inproc_routes(fx.db) do call
                _post(call, "/api/exposures/$e1/assignment/members", Dict(:index_id => 10))
                _post(call, "/api/exposures/$e2/assignment/members", Dict(:index_id => 20))

                @test _post(call, "/api/samples/$(fx.s_id)/call/state", Dict(:state => "bogus")).status == 400
                @test _post(call, "/api/samples/999999/call/state", Dict(:state => "null")).status == 404

                r = _post(call, "/api/samples/$(fx.s_id)/call/state", Dict(:state => "form_factor"))
                @test r.status == 200
                @test _json(r).state == "form_factor"
                @test isempty(_members(fx.db))
                s = _listed(call, fx.exp_id, fx.s_id)
                @test s.assignment_state == "form_factor"
                @test s.phase === nothing
                for e in (e1, e2)
                    @test _json(call("GET", "/api/exposures/$e/assignment")).state == "form_factor"
                end

                _post(call, "/api/exposures/$e2/assignment/members", Dict(:index_id => 20))
                s = _listed(call, fx.exp_id, fx.s_id)
                @test s.assignment_state == "indexed"
                @test s.phase == "HII"
            end
        end
    end
end

@testset "the sample's call state round-trips through the event log" begin
    mktempdir() do dir
        fx = _call_fixture(dir)
        e1, e2 = fx.eids
        _index!(fx.db, 10, e1, "Pn3m", 0.5)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        HimalayaUI.apply_event!(fx.db, req; kind = "assignment_add",
            entity_type = "exposure", entity_id = e1, payload = Dict(:index_id => 10))
        HimalayaUI.apply_event!(fx.db, req; kind = "sample_call_set_state",
            entity_type = "sample", entity_id = fx.s_id, payload = Dict(:state => "null"))
        live = HimalayaUI.sample_call_state(fx.db, fx.s_id)
        @test live == "null"
        @test isempty(_members(fx.db))

        DBInterface.execute(fx.db, "DELETE FROM sample_calls")
        HimalayaUI.rebuild_views_from_log!(fx.db, fx.s_id; entity_type = "sample")
        @test HimalayaUI.sample_call_state(fx.db, fx.s_id) == live
    end
end

@testset "migrate_sample_calls! keeps the call the user made last, and its state" begin
    mktempdir() do dir
        fx = _call_fixture(dir; frames = 3)
        e1, e2, e3 = fx.eids
        _index!(fx.db, 10, e1, "Pn3m", 0.5)
        _index!(fx.db, 20, e2, "Im3m", 0.9)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        for (e, ix) in ((e2, 20), (e1, 10))
            HimalayaUI.apply_event!(fx.db, req; kind = "assignment_add",
                entity_type = "exposure", entity_id = e, payload = Dict(:index_id => ix))
        end
        DBInterface.execute(fx.db, "UPDATE exposures SET status = 'rejected' WHERE id = ?", [e1])
        DBInterface.execute(fx.db,
            "INSERT INTO assignments (exposure_id, state) VALUES (?, 'form_factor')", [e3])
        DBInterface.execute(fx.db, "DELETE FROM sample_calls")
        DBInterface.execute(fx.db, "DELETE FROM schema_migrations WHERE name = ?",
            [HimalayaUI.MIGRATION_SAMPLE_CALLS])

        HimalayaUI.migrate_sample_calls!(fx.db)
        @test _members(fx.db) == [(e1, 10)]
        @test HimalayaUI.sample_call_state(fx.db, fx.s_id) == "indexed"
        @test HimalayaUI._migrated(fx.db, HimalayaUI.MIGRATION_SAMPLE_CALLS)
    end
end

@testset "migrate_sample_calls! without events falls back to the representative frame" begin
    mktempdir() do dir
        fx = _call_fixture(dir)
        e1, e2 = fx.eids
        _index!(fx.db, 10, e1, "Pn3m", 0.5)
        _index!(fx.db, 20, e2, "Im3m", 0.9)
        DBInterface.execute(fx.db, "INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed'), (?, 'indexed')", [e1, e2])
        DBInterface.execute(fx.db, "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, 10), (?, 20)", [e1, e2])
        DBInterface.execute(fx.db, "UPDATE exposures SET selected = 1 WHERE id = ?", [e1])
        DBInterface.execute(fx.db, "DELETE FROM sample_calls")
        DBInterface.execute(fx.db, "DELETE FROM schema_migrations WHERE name = ?",
            [HimalayaUI.MIGRATION_SAMPLE_CALLS])

        HimalayaUI.migrate_sample_calls!(fx.db)
        @test _members(fx.db) == [(e1, 10)]
    end
end
