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

    @testset "a lone live sample leaves when its only partner was merged away" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c, d)) = _seed_bridge_samples(db, 4)
        with_inproc_routes(db) do call
            call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            call("POST", "/api/samples/$c/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(a))
            DBInterface.execute(db, "UPDATE samples SET merged_into_id = ? WHERE id = ?", [d, b])
            call("DELETE", "/api/samples/$c/bridge"; headers=_BRIDGE_HEADERS)
            @test _bridge_rows(db) == Set{Int}()
        end
    end

    @testset "a merged sample cannot be bridged on either side" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c, d)) = _seed_bridge_samples(db, 4)
        with_inproc_routes(db) do call
            call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            DBInterface.execute(db, "UPDATE samples SET merged_into_id = ? WHERE id = ?", [c, b])
            @test call("POST", "/api/samples/$b/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(d)).status == 404
            @test call("POST", "/api/samples/$d/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b)).status == 404
            @test _groups(call, exp_id) == Set()
        end
    end

    @testset "merging a bridged sample hands its bridges to the survivor" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c, d, e)) = _seed_bridge_samples(db, 5)
        merge(loser, survivor) = call -> call("POST", "/api/samples/$loser/merge"; headers=_BRIDGE_HEADERS,
            body=Vector{UInt8}(JSON3.write(Dict(:survivor_id => survivor))))
        with_inproc_routes(db) do call
            call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            @test merge(b, c)(call).status == 200
            @test _groups(call, exp_id) == Set([Set([a, c])])
            call("POST", "/api/samples/$d/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(e))
            @test merge(c, d)(call).status == 200
            @test _groups(call, exp_id) == Set([Set([a, d, e])])
            @test merge(e, a)(call).status == 200
            @test _groups(call, exp_id) == Set([Set([a, d])])
        end
        db2 = open_prepared_clone(mktempdir())
        (exp2, (p, q, r_)) = _seed_bridge_samples(db2, 3)
        with_inproc_routes(db2) do call
            call("POST", "/api/samples/$p/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(q))
            @test merge(r_, p)(call).status == 200
            @test _groups(call, exp2) == Set([Set([p, q])])
            @test _bridge_rows(db2) == Set([p, q])
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

    @testset "a group left with one live sample is not listed" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c)) = _seed_bridge_samples(db, 3)
        with_inproc_routes(db) do call
            call("POST", "/api/samples/$a/bridge"; headers=_BRIDGE_HEADERS, body=_bridge_body(b))
            DBInterface.execute(db, "UPDATE samples SET merged_into_id = ? WHERE id = ?", [c, b])
            @test _groups(call, exp_id) == Set()
        end
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
end
