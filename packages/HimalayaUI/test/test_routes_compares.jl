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

function _seed_run(db, exp_id::Int, load_index::Int, acquisitions)
    load_id = DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO loads (experiment_id, load_index) VALUES (?, ?)", [exp_id, load_index]))
    map(acquisitions) do (kind, scan)
        well = "01" * string('A' + (scan % 8))
        sid = HimalayaUI.create_sample!(db; experiment_id=exp_id, load_id=load_id, slot_index=scan,
            name="lipidA $well ($kind$(lpad(scan, 4, '0')))")
        HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=sid,
            filename="lipidA_$(well)_$(kind)$(lpad(scan, 4, '0'))")
        sid
    end
end

function _seed_autosampler(acquisitions...; mode = "autosampler")
    db = open_prepared_clone(mktempdir())
    exp_id = HimalayaUI.create_experiment!(db; path="/d", data_dir="/d", analysis_dir="/a")
    mode === nothing || DBInterface.execute(db, "UPDATE setups SET mode = ?", [mode])
    runs = [_seed_run(db, exp_id, i, acq) for (i, acq) in enumerate(acquisitions)]
    (db, exp_id, runs)
end

@testset "default background" begin
    run = [("B", 100), ("S", 101), ("S", 102), ("B", 103)]

    @testset "a sample of an autosampler run defaults to the buffer at the start of the run" begin
        (db, _, ((b1, s1, s2, b2),)) = _seed_autosampler(run)
        @test HimalayaUI.default_background(db, s1) == b1
        @test HimalayaUI.default_background(db, s2) == b1
    end

    @testset "the buffer with the lowest scan wins, whatever order the rows were made in" begin
        (db, _, ((s1, b2, b1),)) = _seed_autosampler([("S", 101), ("B", 130), ("B", 100)])
        @test HimalayaUI.default_background(db, s1) == b1
    end

    @testset "a buffer has no default, and nor does a run without a buffer" begin
        (db, _, ((b1, s1, s2, b2), (t1, t2))) = _seed_autosampler(run, [("S", 201), ("S", 202)])
        @test HimalayaUI.default_background(db, b1) === nothing
        @test HimalayaUI.default_background(db, b2) === nothing
        @test HimalayaUI.default_background(db, t1) === nothing
    end

    @testset "a run never borrows the buffer of another run" begin
        (db, _, ((b1, s1, s2, b2), (c1, c2))) = _seed_autosampler(run, [("B", 200), ("S", 201)])
        @test HimalayaUI.default_background(db, c2) == c1
        @test HimalayaUI.default_background(db, s1) == b1
    end

    @testset "capillary and undecided setups have no default" begin
        for mode in ("capillary", nothing)
            (db, _, ((b1, s1, s2, b2),)) = _seed_autosampler(run; mode)
            @test HimalayaUI.default_background(db, s1) === nothing
        end
    end

    @testset "the defaults of a whole experiment come from one pass" begin
        (db, exp_id, ((b1, s1, s2, b2), (c1, c2))) = _seed_autosampler(run, [("B", 200), ("S", 201)])
        @test HimalayaUI.default_backgrounds(db, exp_id) ==
              Dict(b1 => nothing, s1 => b1, s2 => b1, b2 => nothing, c1 => nothing, c2 => c1)
    end
end

_patch(call, sid; user = "alice", kw...) =
    call("PATCH", "/api/samples/$sid/compares"; headers=_compare_headers(user), body=_compare_body(; kw...))

_null_subtraction() = Vector{UInt8}("{\"subtraction\": null}")

_background(body, sid) = only(b for b in body.backgrounds if Int(b.sample_id) == sid)

_stored(db, sid) = [(r.ref_sample_id,) for r in Tables.rowtable(DBInterface.execute(db,
    "SELECT ref_sample_id FROM sample_subtractions WHERE sample_id = ?", [sid]))]

@testset "autosampler backgrounds" begin
    run = [("B", 100), ("S", 101), ("S", 102), ("B", 103)]

    @testset "the experiment reports each autosampler sample's effective background and where it came from" begin
        (db, exp_id, ((b1, s1, s2, b2),)) = _seed_autosampler(run)
        with_inproc_routes(db) do call
            body = _compares(call, exp_id)
            @test _background(body, s1).ref_sample_id == b1
            @test _background(body, s1).source == "default"
            @test _background(body, s1).default_ref_sample_id == b1
            @test _background(body, b1).ref_sample_id === nothing
            @test _background(body, b1).default_ref_sample_id === nothing
            @test isempty(body.subtractions)
            @test _stored(db, s1) == []
        end
    end

    @testset "each entry says whether the sample is a buffer and lists the buffers of its run in scan order" begin
        (db, exp_id, ((b1, s1, s2, b2), (t1,))) = _seed_autosampler(run, [("S", 201)])
        with_inproc_routes(db) do call
            body = _compares(call, exp_id)
            @test _background(body, b2).is_buffer == true
            @test _background(body, s2).is_buffer == false
            @test collect(Int, _background(body, s1).run_buffer_ids) == [b1, b2]
            @test collect(Int, _background(body, b2).run_buffer_ids) == [b1, b2]
            @test collect(Int, _background(body, t1).run_buffer_ids) == Int[]
        end
    end

    @testset "a capillary experiment reports no backgrounds and keeps its subtractions" begin
        (db, exp_id, ((b1, s1, s2, b2),)) = _seed_autosampler(run; mode = "capillary")
        with_inproc_routes(db) do call
            _patch(call, s1; add=[b1], subtraction=b1)
            body = _compares(call, exp_id)
            @test isempty(body.backgrounds)
            @test _sub(body, s1) == b1
        end
    end

    @testset "choosing another sample stores it as the user's choice and uses no comparison slot" begin
        (db, exp_id, ((b1, s1, s2, b2),)) = _seed_autosampler(run)
        with_inproc_routes(db) do call
            r = _patch(call, s1; subtraction=b2)
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test _background(body, s1).ref_sample_id == b2
            @test _background(body, s1).source == "user"
            @test _background(body, s1).default_ref_sample_id == b1
            @test _sub(body, s1) == b2
            @test _refs(body, s1) == Int[]
        end
    end

    @testset "none is stored as its own state, and the default does not come back" begin
        (db, exp_id, ((b1, s1, s2, b2),)) = _seed_autosampler(run)
        with_inproc_routes(db) do call
            _patch(call, s1; subtraction=b2)
            r = call("PATCH", "/api/samples/$s1/compares"; headers=_compare_headers(), body=_null_subtraction())
            @test r.status == 200
            body = JSON3.read(String(r.body))
            @test _background(body, s1).ref_sample_id === nothing
            @test _background(body, s1).source == "user"
            @test _sub(body, s1) === nothing
            @test ismissing(only(_stored(db, s1))[1])
            @test _background(_compares(call, exp_id), s2).ref_sample_id == b1
        end
    end

    @testset "back to the default deletes the stored choice" begin
        (db, exp_id, ((b1, s1, s2, b2),)) = _seed_autosampler(run)
        with_inproc_routes(db) do call
            call("PATCH", "/api/samples/$s1/compares"; headers=_compare_headers(), body=_null_subtraction())
            r = _patch(call, s1; subtraction="default")
            @test r.status == 200
            @test _background(JSON3.read(String(r.body)), s1).ref_sample_id == b1
            @test _background(JSON3.read(String(r.body)), s1).source == "default"
            @test _stored(db, s1) == []
        end
    end

    @testset "adding or removing a comparison leaves the background alone" begin
        (db, exp_id, ((b1, s1, s2, b2),)) = _seed_autosampler(run)
        with_inproc_routes(db) do call
            _patch(call, s1; subtraction=b2)
            _patch(call, s1; add=[s2])
            _patch(call, s1; remove=[s2])
            @test _background(_compares(call, exp_id), s1).ref_sample_id == b2
        end
    end

    @testset "a change to what it already is logs nothing" begin
        (db, exp_id, ((b1, s1, s2, b2),)) = _seed_autosampler(run)
        with_inproc_routes(db) do call
            _patch(call, s1; subtraction="default")
            _patch(call, s1; subtraction=b2)
            _patch(call, s1; subtraction=b2)
        end
        n = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM user_actions WHERE action = 'compares_changed'"))).n)
        @test n == 1
    end

    @testset "the background must be a sample of the same experiment and not the sample itself" begin
        (db, exp_id, ((b1, s1, s2, b2),)) = _seed_autosampler(run)
        (exp2, (x,)) = _seed_compare_samples(db, 1)
        with_inproc_routes(db) do call
            @test _patch(call, s1; subtraction=s1).status == 400
            @test _patch(call, s1; subtraction=x).status == 422
            @test _patch(call, s1; subtraction="buffer").status == 400
        end
    end

    @testset "default is accepted for a capillary sample and clears its subtraction" begin
        (db, exp_id, ((b1, s1, s2, b2),)) = _seed_autosampler(run; mode = "capillary")
        with_inproc_routes(db) do call
            _patch(call, s1; add=[b1], subtraction=b1)
            r = _patch(call, s1; subtraction="default")
            @test r.status == 200
            @test _sub(JSON3.read(String(r.body)), s1) === nothing
        end
    end

    @testset "an event log from before none existed replays the same, and the new states replay too" begin
        (db, exp_id, ((b1, s1, s2, b2),)) = _seed_autosampler(run)
        (a, b) = _seed_compare_samples(db, 2)[2]
        with_inproc_routes(db) do call
            _patch(call, a; add=[b], subtraction=b)
            call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_null_subtraction())
            _patch(call, s1; subtraction=b2)
            call("PATCH", "/api/samples/$s2/compares"; headers=_compare_headers(), body=_null_subtraction())
            _patch(call, s2; subtraction=b2)
            call("PATCH", "/api/samples/$s2/compares"; headers=_compare_headers(), body=_null_subtraction())
        end
        q(sql) = Tables.rowtable(DBInterface.execute(db, sql))
        live = q("SELECT sample_id, ref_sample_id FROM sample_subtractions ORDER BY sample_id")
        @test length(live) == 2
        DBInterface.execute(db, "DELETE FROM sample_subtractions")
        for s in (a, s1, s2)
            HimalayaUI.rebuild_views_from_log!(db, s; entity_type = "sample")
        end
        @test isequal(q("SELECT sample_id, ref_sample_id FROM sample_subtractions ORDER BY sample_id"), live)
    end
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

    @testset "a sample holds at most five comparisons, one per comparison colour" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, rest...)) = _seed_compare_samples(db, 8)
        with_inproc_routes(db) do call
            patch(; kw...) = call("PATCH", "/api/samples/$a/compares"; headers=_compare_headers(), body=_compare_body(; kw...))
            @test patch(add=collect(rest[1:4])).status == 200
            @test patch(add=[rest[5], rest[6]]).status == 422
            @test _refs(_compares(call, exp_id), a) == collect(rest[1:4])
            @test patch(add=[rest[5]]).status == 200
            @test patch(add=[rest[6]]).status == 422
            @test patch(add=[rest[1]]).status == 200
            @test patch(add=[rest[6]], remove=[rest[1]]).status == 200
            @test _refs(_compares(call, exp_id), a) == [rest[2:5]..., rest[6]]
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
        @test cc(Set([1]), 2, [2], Int[], false, nothing) == Dict(:add => [2], :subtraction => nothing)
        @test cc(Set([1]), 2, Int[], [1], false, nothing) == Dict(:remove => [1], :subtraction => nothing)
        @test cc(Set([1]), 2, Int[], Int[], false, nothing) == Dict()
        @test cc(Set([1]), 2, [2], Int[], true, 2) == Dict(:add => [2])
    end

    @testset "compare_changes for an autosampler sample never clears the background with a chip" begin
        cc(args...) = HimalayaUI.compare_changes(args...; managed = true)
        @test cc(Set([1]), 2, Int[], [1], false, nothing) == Dict(:remove => [1])
        @test cc(Set([1]), 2, [3], Int[], false, nothing) == Dict(:add => [3])
        @test cc(Set([1]), nothing, Int[], Int[], true, nothing) == Dict(:subtraction => "none")
        @test cc(Set([1]), :none, Int[], Int[], true, nothing) == Dict()
        @test cc(Set([1]), :none, Int[], Int[], true, "default") == Dict(:subtraction => nothing)
        @test cc(Set([1]), nothing, Int[], Int[], true, "default") == Dict()
        @test cc(Set([1]), 2, Int[], Int[], true, "default") == Dict(:subtraction => nothing)
        @test cc(Set([1]), :none, Int[], Int[], true, 2) == Dict(:subtraction => 2)
        @test cc(Set([1]), 2, Int[], Int[], true, 2) == Dict()
    end

    @testset "compare_changes treats default like null for a capillary sample" begin
        cc = HimalayaUI.compare_changes
        @test cc(Set([1]), 1, Int[], Int[], true, "default") == Dict(:subtraction => nothing)
        @test cc(Set([1]), nothing, Int[], Int[], true, "default") == Dict()
    end
end
