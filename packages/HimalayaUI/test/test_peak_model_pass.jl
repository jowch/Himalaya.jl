using Test, HTTP, JSON3, Random, SQLite, DBInterface, Tables
using HimalayaUI
import HimalayaPeakModel as HPM
isdefined(@__MODULE__, :_pt_db) || include("peak_training_fixture.jl")
isdefined(@__MODULE__, :_pma_db) || include("peak_model_apply_fixture.jl")

_pmp_events(db) = only(Tables.rowtable(DBInterface.execute(db, "SELECT COUNT(*) AS n FROM user_actions"))).n

function _pmp_subscribe()
    sub = (pending = Channel{String}(256),)
    lock(() -> push!(HimalayaUI.SSE_SUBSCRIBERS[], sub), HimalayaUI.SSE_LOCK)
    sub
end

function _pmp_frames(sub)
    out = Any[]
    while isready(sub.pending)
        frame = take!(sub.pending)
        push!(out, JSON3.read(split(split(frame, "data: ")[2], "\n")[1]))
    end
    out
end

function _pmp_with_pass(f)
    _pma_trained()
    HimalayaUI.PEAK_MODEL_PASS.stopped = false
    sub = _pmp_subscribe()
    try
        f(sub)
    finally
        HimalayaUI.stop_peak_model_pass!()
        HimalayaUI.PEAK_MODEL_PREDICT[] = HPM.predict_peaks
        lock(() -> filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[]), HimalayaUI.SSE_LOCK)
        close(sub.pending)
    end
    @test !HimalayaUI.PEAK_MODEL_PASS.running
end

function _pmp_analyzed(dir; kw...)
    b = _pma_db(dir; kw...)
    for e in b.eids
        HimalayaUI.analyze_exposure!(b.db, e, b.adirs[e])
    end
    b
end

@testset "promotion starts a pass that re-detects only unreviewed exposures, quietly" begin
    mktempdir() do dir
        b = _pmp_analyzed(dir; n = 2, experiments = 2)
        reviewed_event, reviewed_flag = b.eids[1], b.eids[3]
        HimalayaUI.apply_event!(b.db, _pma_req(); kind = "peak_added", entity_type = "exposure",
                                entity_id = reviewed_event, payload = Dict(:q => 0.123), defer_broadcast = true)
        _pma_review!(b.db, reviewed_flag)
        pinned = (reviewed_event, reviewed_flag)
        movable = [e for e in b.eids if !(e in pinned)]
        before = Dict(e => _pma_peaks(b.db, e) for e in b.eids)
        events = _pmp_events(b.db)
        versions = Dict(e => _pma_one(b.db, "SELECT analysis_version AS v, trace_hash AS h FROM exposures WHERE id = ?", [e]) for e in b.eids)
        _pmp_with_pass() do sub
            id = _pma_add_model!(b.db)
            @test HimalayaUI.promote_peak_model!(b.db, id)
            @test HimalayaUI.start_peak_model_pass!(b.db) === :started
            HimalayaUI.wait_peak_model_pass()
            for e in movable
                @test _pma_version(b.db, e).v == id
                @test all(p -> p.model_id == id && p.found_by in ("model", "both"), _pma_peaks(b.db, e))
            end
            for e in pinned
                @test isequal(_pma_peaks(b.db, e), before[e])
                @test ismissing(_pma_version(b.db, e).v)
            end
            @test _pmp_events(b.db) == events
            for e in b.eids
                @test _pma_one(b.db, "SELECT analysis_version AS v, trace_hash AS h FROM exposures WHERE id = ?", [e]) == versions[e]
            end
            frames = _pmp_frames(sub)
            @test all(f -> f.kind == "analyze_run" && f.entity_type == "exposure" && f.id == 0, frames)
            @test sort([f.entity_id for f in frames]) == sort(movable)
            @test all(f -> f.payload.findpeaks_skipped == false && haskey(f.post_state, :indices), frames)
            @test first(frames).entity_id == b.eids[4]
            @test all(xid -> !HimalayaUI.ingest_has_work(b.db, xid), b.xids)
        end
        DBInterface.close!(b.db)
    end
end

@testset "a curation that lands during the pass is not overwritten" begin
    mktempdir() do dir
        b = _pmp_analyzed(dir; n = 2)
        target = b.eids[1]
        before = _pma_peaks(b.db, target)
        _pmp_with_pass() do sub
            id = _pma_add_model!(b.db)
            HimalayaUI.promote_peak_model!(b.db, id)
            seen = Ref(false)
            HimalayaUI.PEAK_MODEL_PREDICT[] = function (model, q, I, σ; kw...)
                if !seen[]
                    seen[] = true
                    _pma_review!(b.db, target)
                end
                HPM.predict_peaks(model, q, I, σ; kw...)
            end
            HimalayaUI.start_peak_model_pass!(b.db)
            HimalayaUI.wait_peak_model_pass()
            first_exposure = b.eids[1] == target
            @test isequal(_pma_peaks(b.db, target), before)
            @test ismissing(_pma_version(b.db, target).v)
            @test _pma_version(b.db, b.eids[2]).v == id
            @test [f.entity_id for f in _pmp_frames(sub)] == [b.eids[2]]
        end
        DBInterface.close!(b.db)
    end
end

@testset "a failure on one exposure does not stop the pass" begin
    mktempdir() do dir
        b = _pmp_analyzed(dir; n = 3)
        broken = b.eids[1]
        _pmp_with_pass() do sub
            id = _pma_add_model!(b.db)
            HimalayaUI.promote_peak_model!(b.db, id)
            rm(joinpath(b.adirs[broken], "x1k1_tot.dat"))
            @test_logs (:warn, r"model pass") begin
                HimalayaUI.start_peak_model_pass!(b.db)
                HimalayaUI.wait_peak_model_pass()
            end
            @test ismissing(_pma_version(b.db, broken).v)
            @test ismissing(_pma_one(b.db, "SELECT analysis_error AS e FROM exposures WHERE id = ?", [broken]).e)
            @test all(e -> _pma_version(b.db, e).v == id, b.eids[2:3])
        end
        DBInterface.close!(b.db)
    end
end

@testset "one pass at a time: a start during a pass queues one more sweep" begin
    mktempdir() do dir
        b = _pmp_analyzed(dir; n = 2)
        _pmp_with_pass() do sub
            id = _pma_add_model!(b.db)
            HimalayaUI.promote_peak_model!(b.db, id)
            gate = Channel{Nothing}(1)
            entered = Channel{Nothing}(1)
            HimalayaUI.PEAK_MODEL_PREDICT[] = function (model, q, I, σ; kw...)
                isready(entered) || put!(entered, nothing)
                isopen(gate) && try wait(gate) catch end
                HPM.predict_peaks(model, q, I, σ; kw...)
            end
            @test HimalayaUI.start_peak_model_pass!(b.db) === :started
            take!(entered)
            task = HimalayaUI.PEAK_MODEL_PASS.task
            @test HimalayaUI.start_peak_model_pass!(b.db) === :queued
            @test HimalayaUI.PEAK_MODEL_PASS.task === task
            close(gate)
            HimalayaUI.wait_peak_model_pass()
            @test all(e -> _pma_version(b.db, e).v == id, b.eids)
            @test !HimalayaUI.PEAK_MODEL_PASS.running && !HimalayaUI.PEAK_MODEL_PASS.pending
        end
        DBInterface.close!(b.db)
    end
end

@testset "stopping ends the pass and refuses further starts" begin
    mktempdir() do dir
        b = _pmp_analyzed(dir; n = 4)
        _pmp_with_pass() do sub
            id = _pma_add_model!(b.db)
            HimalayaUI.promote_peak_model!(b.db, id)
            calls = Ref(0)
            HimalayaUI.PEAK_MODEL_PREDICT[] = function (model, q, I, σ; kw...)
                calls[] += 1
                calls[] == 1 && @async HimalayaUI.stop_peak_model_pass!()
                sleep(0.2)
                HPM.predict_peaks(model, q, I, σ; kw...)
            end
            HimalayaUI.start_peak_model_pass!(b.db)
            HimalayaUI.wait_peak_model_pass()
            @test calls[] < 4
            @test HimalayaUI.start_peak_model_pass!(b.db) === :stopped
            @test count(e -> !ismissing(_pma_version(b.db, e).v), b.eids) < 4
        end
        DBInterface.close!(b.db)
    end
end

@testset "curation routes never re-detect because of the model" begin
    mktempdir() do dir
        b = _pmp_analyzed(dir; n = 1)
        e = b.eids[1]
        id = _pma_add_model!(b.db)
        HimalayaUI.promote_peak_model!(b.db, id)
        before = _pma_peaks(b.db, e)
        with_inproc_routes(b.db) do call
            r = call("POST", "/api/exposures/$e/peaks"; headers = ["X-Username" => "alice", "Content-Type" => "application/json",
                     "X-Client-Op-Id" => "op-1"], body = Vector{UInt8}(JSON3.write(Dict(:q => 0.123))))
            @test r.status == 201
        end
        @test isequal(_pma_peaks(b.db, e), before)
        @test ismissing(_pma_version(b.db, e).v)
        DBInterface.close!(b.db)
    end
end

@testset "boot retires a stale model, and resumes the pass for the active one" begin
    mktempdir() do dir
        b = _pmp_analyzed(dir; n = 2)
        _pmp_with_pass() do sub
            stale = _pma_add_model!(b.db; status = "active")
            DBInterface.execute(b.db, "UPDATE peak_models SET feature_version = 'default-1;old' WHERE id = ?", [stale])
            HimalayaUI.boot_peak_model!(b.db)
            @test _pma_status(b.db, stale).status == "retired"
            @test all(e -> ismissing(_pma_version(b.db, e).v), b.eids)
            id = _pma_add_model!(b.db; status = "active")
            HimalayaUI.boot_peak_model!(b.db)
            HimalayaUI.wait_peak_model_pass()
            @test all(e -> _pma_version(b.db, e).v == id, b.eids)
        end
        DBInterface.close!(b.db)
    end
end

@testset "a recorded passed run is promoted and starts the pass" begin
    mktempdir() do dir
        b = _pmp_analyzed(dir; n = 2)
        old = HimalayaUI.PEAK_TRAINING_COMMAND[]
        _pmp_with_pass() do sub
            HimalayaUI.PEAK_TRAINING.stopped = false
            HimalayaUI.PEAK_TRAINING.attempted_watermark = 0
            result = JSON3.write(_pma_result(watermark = 3))
            HimalayaUI.PEAK_TRAINING_COMMAND[] = path -> addenv(`$(Base.julia_cmd()) --startup-file=no -e "println(ENV[\"FAKE_RESULT\"])"`,
                                                                 "FAKE_RESULT" => result)
            try
                for _ in 1:3
                    DBInterface.execute(b.db, "INSERT INTO user_actions (action, entity_type, entity_id, timestamp) VALUES ('peak_added', 'exposure', 1, datetime('now', '-1 hour'))")
                end
                @test HimalayaUI.maybe_start_peak_training!(b.db; min_events = 3, quiet_seconds = 60.0) === :started
                HimalayaUI.wait_peak_training()
                HimalayaUI.wait_peak_model_pass()
                row = only(Tables.rowtable(DBInterface.execute(b.db, "SELECT id, status FROM peak_models")))
                @test row.status == "active"
                @test all(e -> _pma_version(b.db, e).v == row.id, b.eids)
            finally
                HimalayaUI.PEAK_TRAINING_COMMAND[] = old
                HimalayaUI.stop_peak_training!()
            end
        end
        DBInterface.close!(b.db)
    end
end
