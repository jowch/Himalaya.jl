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

function _legacy_event!(db, kind, eid, payload)
    DBInterface.execute(db,
        "INSERT INTO user_actions (action, entity_type, entity_id, payload) VALUES (?, 'exposure', ?, ?)",
        [kind, eid, JSON3.write(payload)])
end

function _legacy!(db; states = Dict(), members = Dict())
    for (e, st) in states
        DBInterface.execute(db, "INSERT INTO assignments (exposure_id, state) VALUES (?, ?)", [e, st])
    end
    for (e, ixs) in members, ix in ixs
        DBInterface.execute(db,
            "INSERT OR IGNORE INTO assignments (exposure_id, state) VALUES (?, 'indexed')", [e])
        DBInterface.execute(db,
            "INSERT INTO assignment_members (exposure_id, index_id) VALUES (?, ?)", [e, ix])
    end
    DBInterface.execute(db, "DELETE FROM sample_calls")
    DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?",
        [HimalayaUI.MIGRATION_SAMPLE_CALLS])
end

_reopen(dir, db) = (close(db); HimalayaUI.open_db(joinpath(dir, "h.db")))

_rebuilt_members(db, eids) = begin
    DBInterface.execute(db, "DELETE FROM assignment_members")
    DBInterface.execute(db, "DELETE FROM sample_calls")
    foreach(e -> HimalayaUI.rebuild_views_from_log!(db, e), eids)
    _members(db)
end

@testset "opening a pre-#350 database migrates it to sample calls" begin
    mktempdir() do dir
        fx = _call_fixture(dir)
        e1, e2 = fx.eids
        _index!(fx.db, 10, e1, "Pn3m", 0.5)
        _index!(fx.db, 20, e2, "Im3m", 0.9)
        @test HimalayaUI._migrated(fx.db, HimalayaUI.MIGRATION_SAMPLE_CALLS)
        _legacy!(fx.db; members = Dict(e1 => [10], e2 => [20]))
        @test !HimalayaUI._migrated(fx.db, HimalayaUI.MIGRATION_SAMPLE_CALLS)

        db = _reopen(dir, fx.db)
        @test HimalayaUI._migrated(db, HimalayaUI.MIGRATION_SAMPLE_CALLS)
        @test _members(db) == [(e2, 20)]
        @test HimalayaUI.sample_call_state(db, fx.s_id) == "indexed"
        close(db)
    end
end

@testset "migrate_sample_calls! keeps the frame edited last, not the lowest id" begin
    mktempdir() do dir
        fx = _call_fixture(dir; frames = 3)
        e1, e2, e3 = fx.eids
        _index!(fx.db, 10, e1, "Pn3m", 0.5)
        _index!(fx.db, 20, e2, "Im3m", 0.9)
        _index!(fx.db, 30, e3, "Ia3d", 0.7)
        _legacy!(fx.db; members = Dict(e1 => [10], e2 => [20], e3 => [30]))
        _legacy_event!(fx.db, "assignment_add", e1, Dict(:index_id => 10))
        _legacy_event!(fx.db, "assignment_add", e3, Dict(:index_id => 30))
        _legacy_event!(fx.db, "assignment_add", e2, Dict(:index_id => 20))

        HimalayaUI.migrate_sample_calls!(fx.db)
        @test _members(fx.db) == [(e2, 20)]
        @test _rebuilt_members(fx.db, fx.eids) == [(e2, 20)]
        @test HimalayaUI.sample_call_state(fx.db, fx.s_id) == "indexed"
    end
end

@testset "migrate_sample_calls! skips a frame whose last edit left it empty" begin
    mktempdir() do dir
        fx = _call_fixture(dir)
        e1, e2 = fx.eids
        _index!(fx.db, 10, e1, "Pn3m", 0.5)
        _index!(fx.db, 20, e2, "Im3m", 0.9)
        _legacy!(fx.db; members = Dict(e1 => [10]), states = Dict(e2 => "indexed"))
        _legacy_event!(fx.db, "assignment_add", e1, Dict(:index_id => 10))
        _legacy_event!(fx.db, "assignment_add", e2, Dict(:index_id => 20))
        _legacy_event!(fx.db, "assignment_remove", e2, Dict(:index_id => 20))

        HimalayaUI.migrate_sample_calls!(fx.db)
        @test _members(fx.db) == [(e1, 10)]
        @test _rebuilt_members(fx.db, fx.eids) == [(e1, 10)]
    end
end

@testset "migrate_sample_calls! keeps a form-factor call made last" begin
    mktempdir() do dir
        fx = _call_fixture(dir)
        e1, e2 = fx.eids
        _index!(fx.db, 10, e1, "Pn3m", 0.5)
        _legacy!(fx.db; members = Dict(e1 => [10]), states = Dict(e2 => "form_factor"))
        _legacy_event!(fx.db, "assignment_add", e1, Dict(:index_id => 10))
        _legacy_event!(fx.db, "assignment_set_state", e2, Dict(:state => "form_factor"))

        HimalayaUI.migrate_sample_calls!(fx.db)
        @test HimalayaUI.sample_call_state(fx.db, fx.s_id) == "form_factor"
        @test isempty(_members(fx.db))
        @test isempty(_rebuilt_members(fx.db, fx.eids))
        @test HimalayaUI.sample_call_state(fx.db, fx.s_id) == "form_factor"
    end
end

@testset "migrate_sample_calls! without events keeps the representative frame, not the lowest id" begin
    mktempdir() do dir
        fx = _call_fixture(dir)
        e1, e2 = fx.eids
        _index!(fx.db, 10, e1, "Pn3m", 0.5)
        _index!(fx.db, 20, e2, "Im3m", 0.9)
        _legacy!(fx.db; members = Dict(e1 => [10], e2 => [20]))

        HimalayaUI.migrate_sample_calls!(fx.db)
        @test _members(fx.db) == [(e2, 20)]
    end
end

@testset "migrate_sample_calls! keeps a legacy selected frame, then the selected column is dropped" begin
    mktempdir() do dir
        fx = _call_fixture(dir)
        e1, e2 = fx.eids
        DBInterface.execute(fx.db, "ALTER TABLE exposures ADD COLUMN selected BOOLEAN DEFAULT FALSE")
        DBInterface.execute(fx.db, "UPDATE exposures SET selected = 1 WHERE id = ?", [e1])
        _index!(fx.db, 10, e1, "Pn3m", 0.5)
        _index!(fx.db, 20, e2, "Im3m", 0.9)
        _legacy!(fx.db; members = Dict(e1 => [10], e2 => [20]))

        HimalayaUI.migrate_sample_calls!(fx.db)
        @test _members(fx.db) == [(e1, 10)]

        HimalayaUI.migrate_drop_exposure_selected!(fx.db)
        @test "selected" ∉ HimalayaUI.cols_of(fx.db, "exposures")
        HimalayaUI.migrate_drop_exposure_selected!(fx.db)
    end
end

@testset "rebuilding one frame from the log reproduces its members and the sample's state" begin
    mktempdir() do dir
        fx = _call_fixture(dir)
        e1, e2 = fx.eids
        _index!(fx.db, 10, e1, "Pn3m", 0.5)
        _index!(fx.db, 20, e2, "Im3m", 0.9)
        req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
        ev(kind, type, id, p) = HimalayaUI.apply_event!(fx.db, req; kind = kind,
            entity_type = type, entity_id = id, payload = p)
        ev("assignment_add", "exposure", e2, Dict(:index_id => 20))
        ev("sample_call_set_state", "sample", fx.s_id, Dict(:state => "form_factor"))
        ev("assignment_add", "exposure", e1, Dict(:index_id => 10))
        live = (_members(fx.db), HimalayaUI.sample_call_state(fx.db, fx.s_id))
        @test live == ([(e1, 10)], "indexed")

        for e in (e2, e1)
            DBInterface.execute(fx.db, "DELETE FROM assignment_members WHERE exposure_id = ?", [e])
            DBInterface.execute(fx.db, "DELETE FROM sample_calls")
            HimalayaUI.rebuild_views_from_log!(fx.db, e)
            @test (_members(fx.db), HimalayaUI.sample_call_state(fx.db, fx.s_id)) == live
        end
    end
end

if isdefined(@__MODULE__, :with_test_server)

    @testset "deleting the experiment leaves no sample call rows behind" begin
        mktempdir() do dir
            fx = _call_fixture(dir)
            with_inproc_routes(fx.db) do call
                _post(call, "/api/samples/$(fx.s_id)/call/state", Dict(:state => "form_factor"))
                @test call("DELETE", "/api/experiments/$(fx.exp_id)";
                    headers = ["X-Username" => "alice"]).status == 200
            end
            @test isempty(Tables.rowtable(DBInterface.execute(fx.db, "SELECT * FROM sample_calls")))
        end
    end

    @testset "the call broadcast carries the sample's state and members" begin
        mktempdir() do dir
            fx = _call_fixture(dir)
            e1, _ = fx.eids
            _index!(fx.db, 10, e1, "Pn3m", 0.5)
            pending = Channel{String}(64)
            sub = (pending = pending,)
            lock(() -> push!(HimalayaUI.SSE_SUBSCRIBERS[], sub), HimalayaUI.SSE_LOCK)
            try
                with_inproc_routes(fx.db) do call
                    _post(call, "/api/exposures/$e1/assignment/members", Dict(:index_id => 10))
                    _post(call, "/api/samples/$(fx.s_id)/call/state", Dict(:state => "indexed"))
                end
            finally
                lock(() -> filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[]), HimalayaUI.SSE_LOCK)
                close(pending)
            end
            call_frames = [JSON3.read(strip(last(split(f, "data: "; limit = 2)))) for f in pending
                           if occursin("\"kind\":\"sample_call_set_state\"", f)]
            @test length(call_frames) == 1
            ps = only(call_frames).post_state
            @test ps.call.state == "indexed"
            @test [m.index_id for m in ps.call.members] == [10]
        end
    end
end
