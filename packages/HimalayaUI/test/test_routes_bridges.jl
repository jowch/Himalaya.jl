using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

if !isdefined(@__MODULE__, :with_test_server)
    include("test_http.jl")
end

function _seed_bridge_samples(db, n::Int)
    exp_id = HimalayaUI.create_experiment!(db; path="/d", data_dir="/d", analysis_dir="/a")
    load_id = DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO loads (experiment_id, load_index) VALUES (?, 1)", [exp_id]))
    ids = [HimalayaUI.create_sample!(db; experiment_id=exp_id, load_id=load_id,
               slot_index=i, name="S$i") for i in 1:n]
    (exp_id, ids)
end

const _BRIDGE_HEADERS = ["Content-Type" => "application/json", "X-Username" => "alice"]

_bridge_body(to) = Vector{UInt8}(JSON3.write(Dict(:to_sample_id => to)))

function _groups(call, exp_id)
    r = call("GET", "/api/experiments/$exp_id/bridges")
    @test r.status == 200
    rows = JSON3.read(String(r.body))
    by = Dict{Int, Set{Int}}()
    for row in rows
        push!(get!(by, Int(row.group_id), Set{Int}()), Int(row.sample_id))
    end
    Set(values(by))
end

_bridge_rows(db) = Set(Int(r.sample_id) for r in Tables.rowtable(DBInterface.execute(db,
    "SELECT sample_id FROM sample_bridges")))

_label_body(label) = Vector{UInt8}(JSON3.write(Dict(:label => label)))

_label_states(db, sid) = Dict(String(r.label) => String(r.state) for r in Tables.rowtable(DBInterface.execute(db,
    "SELECT label, state FROM sample_labels WHERE sample_id = ?", [sid])))

@testset "bridges by setup mode" begin
    function seeded()
        db = open_prepared_clone(mktempdir())
        exp_id = HimalayaUI.create_experiment!(db; path="/d", data_dir="/d", analysis_dir="/a")
        setup(folder, mode) = HimalayaUI.create_setup!(db; experiment_id = exp_id, folder = folder,
            root_dir = "/d/$folder", data_dir = "/d/$folder", analysis_dir = "/a/$folder", mode = mode)
        sample(name, st) = HimalayaUI.create_sample!(db; experiment_id = exp_id, name = name, setup_id = st)
        (; db, exp_id, setup, sample)
    end
    bridge(call, a, b) = call("POST", "/api/samples/$a/bridge"; headers = _BRIDGE_HEADERS, body = _bridge_body(b))
    count(db, table) = Int(only(Tables.rowtable(DBInterface.execute(db, "SELECT COUNT(*) AS c FROM $table"))).c)

    @testset "a capillary and an autosampler sample are not bridged, and nothing is written" begin
        (; db, exp_id, setup, sample) = seeded()
        cap, auto = setup("capillary", "capillary"), setup("solution", "autosampler")
        a, b = sample("lipidA 01B (S0102)", cap), sample("water 01A (B0101)", auto)
        with_inproc_routes(db) do call
            events, cached = count(db, "user_actions"), count(db, "idempotent_responses")
            r = call("POST", "/api/samples/$a/bridge"; headers = vcat(_BRIDGE_HEADERS, ["X-Client-Op-Id" => "op-mixed"]), body = _bridge_body(b))
            @test r.status == 422
            @test JSON3.read(String(r.body)).error == "A capillary sample and an autosampler sample can't be the same sample."
            @test (count(db, "user_actions"), count(db, "idempotent_responses")) == (events, cached)
            @test _groups(call, exp_id) == Set()
        end
    end

    @testset "two autosampler samples are bridged, and so are two capillary ones, and a sample with no setup counts as capillary" begin
        (; db, exp_id, setup, sample) = seeded()
        cap, auto = setup("capillary", "capillary"), setup("solution", "autosampler")
        undecided = setup("other", nothing)
        a1, a2 = sample("water 01A (B0101)", auto), sample("water 02A (B0102)", auto)
        c1, c2, c3 = sample("lipidA 01B (S0102)", cap), sample("lipidA 02B (S0103)", cap), sample("lipidA 03B (S0104)", undecided)
        with_inproc_routes(db) do call
            @test bridge(call, a1, a2).status == 200
            @test bridge(call, c1, c2).status == 200
            @test bridge(call, c2, c3).status == 200
            @test _groups(call, exp_id) == Set([Set([a1, a2]), Set([c1, c2, c3])])
        end
    end

    @testset "joining a group of the other mode is refused" begin
        (; db, exp_id, setup, sample) = seeded()
        cap, auto = setup("capillary", "capillary"), setup("solution", "autosampler")
        a1, a2 = sample("water 01A (B0101)", auto), sample("water 02A (B0102)", auto)
        c1 = sample("lipidA 01B (S0102)", cap)
        with_inproc_routes(db) do call
            @test bridge(call, a1, a2).status == 200
            @test bridge(call, c1, a1).status == 422
            @test bridge(call, a2, c1).status == 422
            @test _groups(call, exp_id) == Set([Set([a1, a2])])
        end
    end
end

@testset "sample bridges" begin
    @testset "bridging two unbridged samples puts them in one group" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c)) = _seed_bridge_samples(db, 3)
        with_inproc_routes(db) do call
            @test _groups(call, exp_id) == Set()
            r = call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            @test r.status == 200
            @test _groups(call, exp_id) == Set([Set([a, b])])
            body = JSON3.read(String(r.body))
            @test Set(Int(x.sample_id) for x in body) == Set([a, b])
        end
    end

    @testset "bridging into an existing group joins it, and two groups combine" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c, d)) = _seed_bridge_samples(db, 4)
        with_inproc_routes(db) do call
            call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            call("POST", "/api/samples/$c/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(a))
            @test _groups(call, exp_id) == Set([Set([a, b, c])])
            call("POST", "/api/samples/$d/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(c))
            @test _groups(call, exp_id) == Set([Set([a, b, c, d])])
        end
        db2 = open_prepared_clone(mktempdir())
        (exp2, (p, q, r_, s)) = _seed_bridge_samples(db2, 4)
        with_inproc_routes(db2) do call
            call("POST", "/api/samples/$p/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(q))
            call("POST", "/api/samples/$r_/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(s))
            @test _groups(call, exp2) == Set([Set([p, q]), Set([r_, s])])
            call("POST", "/api/samples/$q/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(s))
            @test _groups(call, exp2) == Set([Set([p, q, r_, s])])
        end
    end

    @testset "removing a bridge takes only that sample out, and a lone sample leaves too" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c)) = _seed_bridge_samples(db, 3)
        with_inproc_routes(db) do call
            call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            call("POST", "/api/samples/$c/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(a))
            r = call("DELETE", "/api/samples/$c/bridge"; headers=_BRIDGE_HEADERS)
            @test r.status == 200
            @test _groups(call, exp_id) == Set([Set([a, b])])
            @test _bridge_rows(db) == Set([a, b])
            call("DELETE", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS)
            @test _groups(call, exp_id) == Set()
            @test _bridge_rows(db) == Set{Int}()
            r = call("DELETE", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS)
            @test r.status == 200
        end
    end

    @testset "bridging rejects bad input" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b)) = _seed_bridge_samples(db, 2)
        other_exp = HimalayaUI.create_experiment!(db; path="/e", data_dir="/e", analysis_dir="/f")
        z = HimalayaUI.create_sample!(db; experiment_id=other_exp, name="Z")
        with_inproc_routes(db) do call
            @test call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(a)).status == 400
            @test call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS,
                       body=Vector{UInt8}(JSON3.write(Dict(:x => 1)))).status == 400
            @test call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(999_999)).status == 404
            @test call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(z)).status == 422
            @test call("DELETE", "/api/samples/999999/bridge"; headers=_BRIDGE_HEADERS).status == 404
            @test _groups(call, exp_id) == Set()
        end
    end

    @testset "a body that is not JSON is a 400" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b)) = _seed_bridge_samples(db, 2)
        with_inproc_routes(db) do call
            @test call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS,
                       body=Vector{UInt8}("not json")).status == 400
        end
    end

    @testset "deleting the experiment leaves no bridge rows behind" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b)) = _seed_bridge_samples(db, 2)
        with_inproc_routes(db) do call
            call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            @test call("DELETE", "/api/experiments/$exp_id"; headers=_BRIDGE_HEADERS).status == 200
        end
        @test _bridge_rows(db) == Set{Int}()
        @test isempty(Tables.rowtable(DBInterface.execute(db, "PRAGMA foreign_key_check")))
    end

    @testset "bridge events are logged on the sample and rebuild the view" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c)) = _seed_bridge_samples(db, 3)
        with_inproc_routes(db) do call
            call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            call("POST", "/api/samples/$c/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            call("DELETE", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS)
        end
        kinds = [String(r.action) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT action FROM user_actions WHERE entity_type = 'sample' ORDER BY id"))]
        @test kinds == ["sample_bridged", "sample_bridged", "sample_bridged", "sample_unbridged"]
        live = Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id, group_id FROM sample_bridges ORDER BY sample_id"))
        DBInterface.execute(db, "DELETE FROM sample_bridges")
        for s in (a, b, c)
            HimalayaUI.rebuild_views_from_log!(db, s; entity_type = "sample")
        end
        rebuilt = Tables.rowtable(DBInterface.execute(db,
            "SELECT sample_id, group_id FROM sample_bridges ORDER BY sample_id"))
        @test live == rebuilt
        @test Set(Int(r.sample_id) for r in live) == Set([b, c])
    end

    @testset "bridging broadcasts one SSE frame per changed sample" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b)) = _seed_bridge_samples(db, 2)
        sub = (pending = Channel{String}(64),)
        with_inproc_routes(db) do call
            lock(HimalayaUI.SSE_LOCK) do
                push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
            end
            call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
        end
        frames = String[]
        while isready(sub.pending)
            push!(frames, take!(sub.pending))
        end
        bridged = filter(f -> occursin("sample_bridged", f), frames)
        @test length(bridged) == 2
        @test all(f -> occursin("\"experiment_id\":$exp_id", f), bridged)
    end
    @testset "bridged samples share their labels, except flight paths" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c)) = _seed_bridge_samples(db, 3)
        DBInterface.execute(db,
            "INSERT INTO sample_labels (sample_id, label, state) VALUES (?, 'Buffer', 'confirmed'), (?, '37C', 'confirmed'), (?, '1.7 m', 'confirmed'), (?, 'gone', 'removed')",
            [a, b, b, b])
        with_inproc_routes(db) do call
            call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            @test _label_states(db, a) == Dict("Buffer" => "confirmed", "37C" => "confirmed")
            @test _label_states(db, b) == Dict("Buffer" => "confirmed", "37C" => "confirmed",
                                               "1.7 m" => "confirmed", "gone" => "removed")
            @test isempty(HimalayaUI.bridged_label_drift(db))

            call("POST", "/api/samples/$b/labels"; headers=_BRIDGE_HEADERS, body=_label_body("Calibrant"))
            @test _label_states(db, a)["Calibrant"] == "confirmed"
            call("DELETE", "/api/samples/$a/labels?label=Buffer"; headers=_BRIDGE_HEADERS)
            @test _label_states(db, b)["Buffer"] == "removed"
            call("POST", "/api/samples/$a/labels/suggested"; headers=_BRIDGE_HEADERS, body=_label_body("Buffer"))
            @test _label_states(db, b)["Buffer"] == "suggested"
            call("POST", "/api/samples/$a/labels"; headers=_BRIDGE_HEADERS, body=_label_body("0.5 m"))
            @test !haskey(_label_states(db, b), "0.5 m")
            @test isempty(HimalayaUI.bridged_label_drift(db))

            call("POST", "/api/samples/$c/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            @test _label_states(db, c) == Dict("37C" => "confirmed", "Calibrant" => "confirmed")
            @test isempty(HimalayaUI.bridged_label_drift(db))
        end
    end

    @testset "leaving a bridge keeps a copy of the labels, and later edits stay put" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b)) = _seed_bridge_samples(db, 2)
        with_inproc_routes(db) do call
            call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            call("POST", "/api/samples/$a/labels"; headers=_BRIDGE_HEADERS, body=_label_body("37C"))
            call("DELETE", "/api/samples/$b/bridge"; headers=_BRIDGE_HEADERS)
            @test _label_states(db, b) == Dict("37C" => "confirmed")
            call("POST", "/api/samples/$a/labels"; headers=_BRIDGE_HEADERS, body=_label_body("Buffer"))
            @test !haskey(_label_states(db, b), "Buffer")
        end
    end

    @testset "shared labels survive a rebuild from the log" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b)) = _seed_bridge_samples(db, 2)
        DBInterface.execute(db, "INSERT INTO sample_labels (sample_id, label, state) VALUES (?, 'Buffer', 'confirmed')", [a])
        with_inproc_routes(db) do call
            call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            call("POST", "/api/samples/$b/labels"; headers=_BRIDGE_HEADERS, body=_label_body("37C"))
        end
        live = _label_states(db, b)
        DBInterface.execute(db, "DELETE FROM sample_labels WHERE sample_id = ?", [b])
        HimalayaUI.rebuild_views_from_log!(db, b; entity_type = "sample")
        @test _label_states(db, b) == live
    end

    @testset "the export names the other bridged sample and its setup, never the row itself" begin
        db = open_prepared_clone(mktempdir())
        exp_id = HimalayaUI.create_experiment!(db; path="/d", data_dir="/d", analysis_dir="/a")
        s1 = HimalayaUI.create_setup!(db; experiment_id=exp_id, data_dir="/d/setup1/data", analysis_dir="/a1", folder="setup1")
        s2 = HimalayaUI.create_setup!(db; experiment_id=exp_id, data_dir="/d/setup2/data", analysis_dir="/a2", folder="setup2", display_name="long pipe")
        a = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="DOPE_H2O", setup_id=s1)
        b = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="DOPE_H2O", setup_id=s2)
        c = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="lone", setup_id=s1)
        with_inproc_routes(db) do call
            call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            body = JSON3.read(String(call("GET", "/api/experiments/$exp_id/export?format=json").body))
            rows = Dict(String(s.setup) * "/" * String(s.name) => s for s in body)
            @test [(String(x.name), String(x.setup)) for x in rows["setup1/DOPE_H2O"].bridged_with] == [("DOPE_H2O", "long pipe")]
            @test [(String(x.name), String(x.setup)) for x in rows["long pipe/DOPE_H2O"].bridged_with] == [("DOPE_H2O", "setup1")]
            @test isempty(rows["setup1/lone"].bridged_with)
            csv = String(call("GET", "/api/experiments/$exp_id/export?format=csv").body)
            header = split(first(split(csv, '\n')), ',')
            @test header[end-2:end] == ["setup", "bridged_with", "analyzed"]
            lines = [split(l, ',') for l in split(strip(csv), '\n')[2:end]]
            by = Dict((l[1], l[end-2]) => l[end-1] for l in lines)
            @test by[("DOPE_H2O", "setup1")] == "DOPE_H2O (long pipe)"
            @test by[("DOPE_H2O", "long pipe")] == "DOPE_H2O (setup1)"
            @test by[("lone", "setup1")] == ""
        end
    end

    @testset "a one-time repair shares labels across bridges made before labels were shared, through the log" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c)) = _seed_bridge_samples(db, 3)
        DBInterface.execute(db, "INSERT INTO sample_bridges (sample_id, group_id) VALUES (?, 9), (?, 9)", [a, b])
        DBInterface.execute(db,
            "INSERT INTO sample_labels (sample_id, label, state) VALUES (?, 'Buffer', 'confirmed'), (?, '37C', 'confirmed'), (?, 'buffer', 'removed'), (?, '1.7 m', 'confirmed'), (?, 'other', 'confirmed')",
            [a, b, b, b, c])
        DBInterface.execute(db, "DELETE FROM schema_migrations WHERE name = ?", [HimalayaUI.MIGRATION_BRIDGED_LABELS])
        @test HimalayaUI.bridged_label_drift(db) == [9]
        HimalayaUI.migrate_share_bridged_labels!(db)
        @test isempty(HimalayaUI.bridged_label_drift(db))
        @test _label_states(db, a) == Dict("Buffer" => "confirmed", "37C" => "confirmed")
        @test _label_states(db, b)["buffer"] == "confirmed"
        @test !haskey(_label_states(db, a), "1.7 m")
        @test _label_states(db, c) == Dict("other" => "confirmed")
        logged = Tables.rowtable(DBInterface.execute(db,
            "SELECT entity_id, payload, user_id FROM user_actions WHERE action = 'label_added' ORDER BY id"))
        @test Set(Int(r.entity_id) for r in logged) == Set([a, b])
        @test all(r -> r.user_id === missing, logged)
        n = length(logged)
        HimalayaUI.migrate_share_bridged_labels!(db)
        @test length(Tables.rowtable(DBInterface.execute(db, "SELECT id FROM user_actions WHERE action = 'label_added'"))) == n
        DBInterface.execute(db, "DELETE FROM sample_labels WHERE sample_id = ?", [a])
        HimalayaUI.rebuild_views_from_log!(db, a; entity_type = "sample")
        @test _label_states(db, a) == Dict("37C" => "confirmed")
    end
end
