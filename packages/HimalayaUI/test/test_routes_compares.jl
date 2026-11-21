using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

if !isdefined(@__MODULE__, :with_test_server)
    include("test_http.jl")
end

function _seed_compare_samples(db, n::Int)
    exp_id = HimalayaUI.create_experiment!(db; path="/d", data_dir="/d", analysis_dir="/a")
    load_id = DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO loads (experiment_id, load_index) VALUES (?, 1)", [exp_id]))
    ids = [HimalayaUI.create_sample!(db; experiment_id=exp_id, load_id=load_id,
               slot_index=i, name="S$i") for i in 1:n]
    (exp_id, ids)
end

_compare_headers(user = "alice") = ["Content-Type" => "application/json", "X-Username" => user]
_compare_body(; kw...) = Vector{UInt8}(JSON3.write(Dict(kw...)))

function _compares(call, exp_id)
    r = call("GET", "/api/experiments/$exp_id/compares")
    @test r.status == 200
    JSON3.read(String(r.body))
end

_refs(body, sid) = [Int(c.ref_sample_id) for c in body.compares if Int(c.sample_id) == sid]

function _sub(body, sid)
    rows = [Int(t.ref_sample_id) for t in body.subtractions if Int(t.sample_id) == sid]
    isempty(rows) ? nothing : only(rows)
end

@testset "sample compares" begin
    @testset "adding and removing comparisons keeps them per sample, in the order added" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c, d)) = _seed_compare_samples(db, 4)
        with_inproc_routes(db) do call
            r = call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(add=[c, b]))
            @test r.status == 200
            @test _refs(JSON3.read(String(r.body)), a) == [c, b]
            call("PATCH", "/api/samples/$d/compares"; headers=_compare_headers(), body=_compare_body(add=[b]))
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(remove=[c]))
            body = _compares(call, exp_id)
            @test _refs(body, a) == [b]
            @test _refs(body, d) == [b]
            @test _refs(body, b) == Int[]
        end
    end

    @testset "each comparison remembers who added it" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c)) = _seed_compare_samples(db, 3)
        with_inproc_routes(db) do call
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers("alice"), body=_compare_body(add=[b]))
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers("bob"), body=_compare_body(add=[c, b]))
            by = Dict(Int(x.ref_sample_id) => x.added_by for x in _compares(call, exp_id).compares)
            @test by == Dict(b => "alice", c => "bob")
        end
    end

    @testset "the subtraction is saved per sample, cleared with null, and cleared when its comparison goes" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c)) = _seed_compare_samples(db, 3)
        with_inproc_routes(db) do call
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(add=[b, c], subtraction=b))
            @test _sub(_compares(call, exp_id), a) == b
            @test _sub(_compares(call, exp_id), b) === nothing
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(subtraction=c))
            @test _sub(_compares(call, exp_id), a) == c
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(remove=[b]))
            @test _sub(_compares(call, exp_id), a) == c
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(remove=[c]))
            @test _sub(_compares(call, exp_id), a) === nothing
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(add=[b], subtraction=b))
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(),
                 body=Vector{UInt8}("{\"subtraction\": null}"))
            body = _compares(call, exp_id)
            @test _sub(body, a) === nothing
            @test _refs(body, a) == [b]
        end
    end

    @testset "a change that changes nothing logs nothing" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c)) = _seed_compare_samples(db, 3)
        with_inproc_routes(db) do call
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(add=[b], subtraction=b))
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(add=[b], subtraction=b))
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(remove=[c]))
        end
        n = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM user_actions WHERE action = 'compares_changed'"))).n)
        @test n == 1
    end

    @testset "bad requests are refused" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b)) = _seed_compare_samples(db, 2)
        (exp2, (x,)) = _seed_compare_samples(db, 1)
        with_inproc_routes(db) do call
            patch(body) = call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=body)
            @test patch(Vector{UInt8}("nope")).status == 400
            @test patch(_compare_body(add="b")).status == 400
            @test patch(_compare_body(add=[a])).status == 400
            @test patch(_compare_body(subtraction=a)).status == 400
            @test patch(_compare_body(subtraction="b")).status == 400
            @test patch(_compare_body(add=[999999])).status == 404
            @test patch(_compare_body(add=[x])).status == 422
            @test call("PATCH", "/api/samples/999999/compares"; headers=_compare_headers(), body=_compare_body(add=[b])).status == 404
            @test _compares(call, exp_id).compares == []
        end
    end

    @testset "merged samples drop out of the list" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c)) = _seed_compare_samples(db, 3)
        with_inproc_routes(db) do call
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(add=[b, c], subtraction=b))
            DBInterface.execute(db, "UPDATE samples SET merged_into_id = ? WHERE id = ?", [c, b])
            body = _compares(call, exp_id)
            @test _refs(body, a) == [c]
            @test _sub(body, a) === nothing
        end
    end

    @testset "compare events are logged on the sample and rebuild the view" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b, c)) = _seed_compare_samples(db, 3)
        with_inproc_routes(db) do call
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(add=[b, c], subtraction=c))
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers("bob"), body=_compare_body(remove=[b]))
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers("bob"), body=_compare_body(add=[b]))
            call("PATCH", "/api/samples/$b/compares"; headers=_compare_headers(), body=_compare_body(add=[a]))
        end
        q(sql) = Tables.rowtable(DBInterface.execute(db, sql))
        live = (q("SELECT sample_id, ref_sample_id, added_by FROM sample_compares ORDER BY sample_id, id"),
                q("SELECT sample_id, ref_sample_id FROM sample_subtractions ORDER BY sample_id"))
        DBInterface.execute(db, "DELETE FROM sample_compares")
        DBInterface.execute(db, "DELETE FROM sample_subtractions")
        for s in (a, b, c)
            HimalayaUI.rebuild_views_from_log!(db, s; entity_type = "sample")
        end
        rebuilt = (q("SELECT sample_id, ref_sample_id, added_by FROM sample_compares ORDER BY sample_id, id"),
                   q("SELECT sample_id, ref_sample_id FROM sample_subtractions ORDER BY sample_id"))
        @test live == rebuilt
        @test [Int(r.ref_sample_id) for r in live[1] if Int(r.sample_id) == a] == [c, b]
        @test length(live[2]) == 1
    end

    @testset "a change broadcasts one SSE frame naming the experiment" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b)) = _seed_compare_samples(db, 2)
        sub = (pending = Channel{String}(64),)
        with_inproc_routes(db) do call
            lock(HimalayaUI.SSE_LOCK) do
                push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
            end
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(add=[b]))
        end
        frames = String[]
        while isready(sub.pending)
            push!(frames, take!(sub.pending))
        end
        changed = filter(f -> occursin("compares_changed", f), frames)
        @test length(changed) == 1
        @test occursin("\"experiment_id\":$exp_id", only(changed))
    end

    @testset "compare_changes keeps only real changes" begin
        cc = HimalayaUI.compare_changes
        @test cc(Set([1]), nothing, [1, 2, 2], Int[], false, nothing) == Dict(:add => [2])
        @test cc(Set([1, 2]), 2, Int[], [2, 3], false, nothing) == Dict(:remove => [2], :subtraction => nothing)
        @test cc(Set([1, 2]), 2, Int[], [2], true, 1) == Dict(:remove => [2], :subtraction => 1)
        @test cc(Set([1]), 1, Int[], Int[], true, 1) == Dict()
        @test cc(Set{Int}(), nothing, [3], [3], false, nothing) == Dict()
    end
end
