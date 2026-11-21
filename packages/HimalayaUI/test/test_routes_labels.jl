using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

if !isdefined(@__MODULE__, :with_test_server)
    include("test_http.jl")
end

function _seed_label_samples(db, n::Int)
    exp_id = HimalayaUI.create_experiment!(db; path="/d", data_dir="/d", analysis_dir="/a")
    load_id = DBInterface.lastrowid(DBInterface.execute(db,
        "INSERT INTO loads (experiment_id, load_index) VALUES (?, 1)", [exp_id]))
    ids = [HimalayaUI.create_sample!(db; experiment_id=exp_id, load_id=load_id,
               slot_index=i, name="S$i") for i in 1:n]
    (exp_id, ids)
end

const _LABEL_HEADERS = ["Content-Type" => "application/json", "X-Username" => "alice"]

_label_body(label) = Vector{UInt8}(JSON3.write(Dict(:label => label)))

_add_label(call, sid, label) =
    call("POST", "/api/samples/$sid/labels"; headers=_LABEL_HEADERS, body=_label_body(label))

_remove_label(call, sid, label) =
    call("DELETE", "/api/samples/$sid/labels?label=" * HTTP.escapeuri(label); headers=_LABEL_HEADERS)

_restore_suggestion(call, sid, label) =
    call("POST", "/api/samples/$sid/labels/suggested"; headers=_LABEL_HEADERS, body=_label_body(label))

_suggest!(db, sid, label) = DBInterface.execute(db,
    "INSERT INTO sample_labels (sample_id, label, state, source) VALUES (?, ?, 'suggested', 'ingest')", [sid, label])

_label_events(db, kind) = length(Tables.rowtable(DBInterface.execute(db,
    "SELECT id FROM user_actions WHERE action = ?", [kind])))

function _listed(call, exp_id)
    r = call("GET", "/api/experiments/$exp_id/labels")
    @test r.status == 200
    Set((Int(x.sample_id), String(x.label), String(x.state)) for x in JSON3.read(String(r.body)))
end

_label_rows(db) = Set((Int(r.sample_id), String(r.label), String(r.state)) for r in
    Tables.rowtable(DBInterface.execute(db, "SELECT sample_id, label, state FROM sample_labels")))

@testset "sample labels" begin
    @testset "adding a label lists it as confirmed for that sample" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b)) = _seed_label_samples(db, 2)
        with_inproc_routes(db) do call
            @test _listed(call, exp_id) == Set()
            r = _add_label(call, a, "37C")
            @test r.status == 200
            @test Set((Int(x.sample_id), String(x.label)) for x in JSON3.read(String(r.body))) ==
                  Set([(a, "37C")])
            _add_label(call, b, "buffer")
            @test _listed(call, exp_id) == Set([(a, "37C", "confirmed"), (b, "buffer", "confirmed")])
        end
    end

    @testset "labels are trimmed and unique per sample ignoring case" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a,)) = _seed_label_samples(db, 1)
        with_inproc_routes(db) do call
            _add_label(call, a, "  Buffer ")
            _add_label(call, a, "buffer")
            @test _listed(call, exp_id) == Set([(a, "Buffer", "confirmed")])
        end
    end

    @testset "removing a label hides it and keeps a removed record" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a,)) = _seed_label_samples(db, 1)
        with_inproc_routes(db) do call
            _add_label(call, a, "retake")
            r = _remove_label(call, a, "RETAKE")
            @test r.status == 200
            @test isempty(JSON3.read(String(r.body)))
            @test _listed(call, exp_id) == Set()
        end
        @test _label_rows(db) == Set([(a, "retake", "removed")])
    end

    @testset "adding a removed label brings it back as confirmed" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a,)) = _seed_label_samples(db, 1)
        with_inproc_routes(db) do call
            _add_label(call, a, "retake")
            _remove_label(call, a, "retake")
            _add_label(call, a, "retake")
            @test _listed(call, exp_id) == Set([(a, "retake", "confirmed")])
        end
    end

    @testset "adding a suggested label confirms it" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a,)) = _seed_label_samples(db, 1)
        DBInterface.execute(db,
            "INSERT INTO sample_labels (sample_id, label, state, source) VALUES (?, 'buffer', 'suggested', 'ingest')", [a])
        with_inproc_routes(db) do call
            @test _listed(call, exp_id) == Set([(a, "buffer", "suggested")])
            _add_label(call, a, "buffer")
            @test _listed(call, exp_id) == Set([(a, "buffer", "confirmed")])
        end
        @test only(Tables.rowtable(DBInterface.execute(db, "SELECT source FROM sample_labels WHERE sample_id = ?", [a]))).source == "user"
    end

    @testset "removing a label the sample does not carry is a no-op without an event" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a,)) = _seed_label_samples(db, 1)
        with_inproc_routes(db) do call
            @test _remove_label(call, a, "nothing").status == 200
        end
        @test isempty(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM user_actions WHERE action = 'label_removed'")))
    end

    @testset "invalid requests are refused" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a,)) = _seed_label_samples(db, 1)
        with_inproc_routes(db) do call
            @test _add_label(call, a, "   ").status == 400
            @test _add_label(call, a, repeat("x", 65)).status == 400
            @test _add_label(call, a, "two\nlines").status == 400
            @test call("POST", "/api/samples/$a/labels"; headers=_LABEL_HEADERS,
                       body=Vector{UInt8}("{}")).status == 400
            @test call("DELETE", "/api/samples/$a/labels"; headers=_LABEL_HEADERS).status == 400
            @test _add_label(call, 999_999, "x").status == 404
            @test _remove_label(call, 999_999, "x").status == 404
        end
        @test _label_rows(db) == Set()
    end

    @testset "label events are logged on the sample and rebuild the view" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b)) = _seed_label_samples(db, 2)
        with_inproc_routes(db) do call
            _add_label(call, a, "37C")
            _add_label(call, a, "retake")
            _add_label(call, b, "buffer")
            _remove_label(call, a, "retake")
        end
        kinds = [String(r.action) for r in Tables.rowtable(DBInterface.execute(db,
            "SELECT action FROM user_actions WHERE entity_type = 'sample' ORDER BY id"))]
        @test kinds == ["label_added", "label_added", "label_added", "label_removed"]
        live = _label_rows(db)
        DBInterface.execute(db, "DELETE FROM sample_labels")
        for s in (a, b)
            HimalayaUI.rebuild_views_from_log!(db, s; entity_type = "sample")
        end
        @test _label_rows(db) == live
        @test live == Set([(a, "37C", "confirmed"), (a, "retake", "removed"), (b, "buffer", "confirmed")])
    end

    @testset "a removed suggestion is rebuilt from the log as removed" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a,)) = _seed_label_samples(db, 1)
        DBInterface.execute(db,
            "INSERT INTO sample_labels (sample_id, label, state, source) VALUES (?, 'buffer', 'suggested', 'ingest')", [a])
        with_inproc_routes(db) do call
            _remove_label(call, a, "buffer")
        end
        DBInterface.execute(db, "DELETE FROM sample_labels")
        HimalayaUI.rebuild_views_from_log!(db, a; entity_type = "sample")
        @test _label_rows(db) == Set([(a, "buffer", "removed")])
    end

    @testset "a confirmed or dismissed suggestion goes back to suggested, from the ingest" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a, b)) = _seed_label_samples(db, 2)
        _suggest!(db, a, "buffer")
        _suggest!(db, b, "calibrant")
        with_inproc_routes(db) do call
            _add_label(call, a, "Buffer")
            _remove_label(call, b, "calibrant")
            r = _restore_suggestion(call, a, "BUFFER")
            @test r.status == 200
            @test (a, "buffer", "suggested") in Set((Int(x.sample_id), String(x.label), String(x.state)) for x in JSON3.read(String(r.body)))
            _restore_suggestion(call, b, "calibrant")
            @test _listed(call, exp_id) == Set([(a, "buffer", "suggested"), (b, "calibrant", "suggested")])
        end
        @test Set(String(r.source) for r in Tables.rowtable(DBInterface.execute(db, "SELECT source FROM sample_labels"))) == Set(["ingest"])
        @test _label_events(db, "label_suggestion_restored") == 2
    end

    @testset "restoring a suggestion is a no-op without an event when the label is absent or already suggested" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a,)) = _seed_label_samples(db, 1)
        _suggest!(db, a, "buffer")
        with_inproc_routes(db) do call
            @test _restore_suggestion(call, a, "buffer").status == 200
            @test _restore_suggestion(call, a, "retake").status == 200
            @test _restore_suggestion(call, a, "  ").status == 400
            @test _restore_suggestion(call, 999_999, "buffer").status == 404
        end
        @test _label_events(db, "label_suggestion_restored") == 0
        @test _label_rows(db) == Set([(a, "buffer", "suggested")])
    end

    @testset "a restored suggestion rebuilds from the log and broadcasts with the experiment" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a,)) = _seed_label_samples(db, 1)
        _suggest!(db, a, "buffer")
        sub = (pending = Channel{String}(64),)
        with_inproc_routes(db) do call
            _add_label(call, a, "buffer")
            lock(HimalayaUI.SSE_LOCK) do
                push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
            end
            _restore_suggestion(call, a, "buffer")
        end
        frames = String[]
        while isready(sub.pending)
            push!(frames, take!(sub.pending))
        end
        restored = filter(f -> occursin("label_suggestion_restored", f), frames)
        @test length(restored) == 1
        @test occursin("\"experiment_id\":$exp_id", only(restored))
        DBInterface.execute(db, "DELETE FROM sample_labels")
        HimalayaUI.rebuild_views_from_log!(db, a; entity_type = "sample")
        @test only(Tables.rowtable(DBInterface.execute(db, "SELECT label, state, source FROM sample_labels"))) ==
              (label = "buffer", state = "suggested", source = "ingest")
    end

    @testset "adding a label broadcasts one SSE frame carrying the experiment" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a,)) = _seed_label_samples(db, 1)
        sub = (pending = Channel{String}(64),)
        with_inproc_routes(db) do call
            lock(HimalayaUI.SSE_LOCK) do
                push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
            end
            _add_label(call, a, "37C")
        end
        frames = String[]
        while isready(sub.pending)
            push!(frames, take!(sub.pending))
        end
        added = filter(f -> occursin("label_added", f), frames)
        @test length(added) == 1
        @test occursin("\"experiment_id\":$exp_id", only(added))
        @test occursin("\"label\":\"37C\"", only(added))
    end

    @testset "deleting the experiment leaves no label rows behind" begin
        db = open_prepared_clone(mktempdir())
        (exp_id, (a,)) = _seed_label_samples(db, 1)
        with_inproc_routes(db) do call
            _add_label(call, a, "37C")
            @test call("DELETE", "/api/experiments/$exp_id"; headers=_LABEL_HEADERS).status == 200
        end
        @test _label_rows(db) == Set()
        @test isempty(Tables.rowtable(DBInterface.execute(db, "PRAGMA foreign_key_check")))
    end
end
