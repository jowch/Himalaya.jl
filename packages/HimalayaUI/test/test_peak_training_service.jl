using Test, JSON3, SQLite, DBInterface, Tables
using HimalayaUI
import HimalayaPeakModel as HPM
isdefined(@__MODULE__, :_pt_db) || include("peak_training_fixture.jl")

_pts_rows(db) = Tables.rowtable(DBInterface.execute(db, "SELECT * FROM peak_models ORDER BY id"))

function _pts_events!(db, n; age = "-1 hour")
    for _ in 1:n
        DBInterface.execute(db,
            "INSERT INTO user_actions (action, entity_type, entity_id, timestamp) VALUES ('peak_added', 'exposure', 1, datetime('now', ?))", [age])
    end
end

_pts_result(; status = "passed", watermark = 0) = JSON3.write(Dict(
    "status" => status, "watermark" => watermark, "counts" => Dict("curations" => 1),
    "recipe" => Dict("version" => 1), "recipe_id" => "r1", "feature_version" => "fv",
    "core_version" => "0.0.0", "threshold" => 0.5, "nms" => 8,
    "metrics" => Dict("gate" => Dict("passed" => status == "passed", "reasons" => String[])),
    "fitted" => status == "passed" ? Dict("kind" => "fake") : nothing))

function _pts_fake(code; env = Pair{String,String}[])
    path -> addenv(`$(Base.julia_cmd()) --startup-file=no -e $code`, env...)
end

function _pts_with_service(f, dir; command)
    db = open_prepared_clone(dir)
    old = HimalayaUI.PEAK_TRAINING_COMMAND[]
    HimalayaUI.PEAK_TRAINING_COMMAND[] = command
    HimalayaUI.PEAK_TRAINING.stopped = false
    HimalayaUI.PEAK_TRAINING.attempted_watermark = 0
    try
        f(db)
    finally
        HimalayaUI.stop_peak_training!()
        HimalayaUI.PEAK_TRAINING_COMMAND[] = old
        DBInterface.close!(db)
    end
end

const _PTS_PRINT = "println(ENV[\"FAKE_RESULT\"])"

@testset "a due run records the child's result, and a quiet database starts nothing" begin
    mktempdir() do dir
        _pts_with_service(dir; command = _pts_fake(_PTS_PRINT; env = ["FAKE_RESULT" => _pts_result(watermark = 4)])) do db
            kw = (min_events = 3, quiet_seconds = 60.0)
            _pts_events!(db, 2)
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :not_due
            @test !HimalayaUI.peak_training_running()
            _pts_events!(db, 1)
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :started
            HimalayaUI.wait_peak_training()
            row = only(_pts_rows(db))
            @test row.status == "passed" && row.label_watermark == 4 && row.recipe_id == "r1"
            @test JSON3.read(row.fitted_json).kind == "fake"
            @test !HimalayaUI.peak_training_running()
            @test HimalayaUI.maybe_start_peak_training!(db; min_events = 1, quiet_seconds = 60.0) === :not_due
            @test length(_pts_rows(db)) == 1
        end
    end
end

@testset "an in-memory database has nowhere for the child to read" begin
    db = SQLite.DB()
    @test HimalayaUI.maybe_start_peak_training!(db; min_events = 0, quiet_seconds = 0.0) === :unavailable
    DBInterface.close!(db)
end

@testset "single flight: a due check during a run starts nothing, and the next run follows the new labels" begin
    mktempdir() do dir
        slow = "sleep(2); println(ENV[\"FAKE_RESULT\"])"
        _pts_with_service(dir; command = _pts_fake(slow; env = ["FAKE_RESULT" => _pts_result(status = "rejected", watermark = 3)])) do db
            kw = (min_events = 3, quiet_seconds = 0.0)
            _pts_events!(db, 3)
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :started
            task = HimalayaUI.PEAK_TRAINING.task
            @test HimalayaUI.peak_training_running()
            _pts_events!(db, 3)
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :busy
            @test HimalayaUI.PEAK_TRAINING.task === task
            HimalayaUI.wait_peak_training()
            row = only(_pts_rows(db))
            @test row.status == "rejected" && ismissing(row.fitted_json)
            @test HimalayaUI.peak_training_due(db; kw...)
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :started
            HimalayaUI.wait_peak_training()
            @test length(_pts_rows(db)) == 2
        end
    end
end

@testset "a crashing child leaves a failed row with the end of its stderr, and is not retried" begin
    mktempdir() do dir
        crash = "print(stderr, \"x\" ^ 5000, \"the real reason\"); exit(3)"
        _pts_with_service(dir; command = _pts_fake(crash)) do db
            kw = (min_events = 2, quiet_seconds = 0.0)
            _pts_events!(db, 2)
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :started
            HimalayaUI.wait_peak_training()
            row = only(_pts_rows(db))
            @test row.status == "failed" && row.label_watermark == 2 && ismissing(row.fitted_json)
            err = JSON3.read(row.metrics_json).error
            @test occursin("the real reason", err) && occursin("code 3", err) && length(err) <= 2000
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :not_due
        end
    end
end

@testset "output that is not a result is a failed row" begin
    mktempdir() do dir
        _pts_with_service(dir; command = _pts_fake("println(\"not json\")")) do db
            _pts_events!(db, 1)
            @test HimalayaUI.maybe_start_peak_training!(db; min_events = 1, quiet_seconds = 0.0) === :started
            HimalayaUI.wait_peak_training()
            row = only(_pts_rows(db))
            @test row.status == "failed" && occursin("no usable result", JSON3.read(row.metrics_json).error)
        end
    end
end

@testset "stopping kills the child and records nothing" begin
    mktempdir() do dir
        _pts_with_service(dir; command = _pts_fake("sleep(60)")) do db
            _pts_events!(db, 1)
            @test HimalayaUI.maybe_start_peak_training!(db; min_events = 1, quiet_seconds = 0.0) === :started
            sleep(0.5)
            proc = HimalayaUI.PEAK_TRAINING.process
            @test proc !== nothing && Base.process_running(proc)
            HimalayaUI.stop_peak_training!()
            @test !Base.process_running(proc) && !HimalayaUI.peak_training_running()
            @test isempty(_pts_rows(db))
        end
    end
end

@testset "the timer checks dueness and stops when told to" begin
    mktempdir() do dir
        _pts_with_service(dir; command = _pts_fake(_PTS_PRINT; env = ["FAKE_RESULT" => _pts_result(watermark = 1)])) do db
            _pts_events!(db, 1)
            HimalayaUI.start_peak_training_timer!(db; interval = 0.2, initial_delay = 0.05,
                                                  min_events = 1, quiet_seconds = 0.0)
            deadline = time() + 60
            while isempty(_pts_rows(db)) && time() < deadline; sleep(0.1); end
            HimalayaUI.wait_peak_training()
            @test length(_pts_rows(db)) == 1
            sleep(0.8)
            @test length(_pts_rows(db)) == 1
            HimalayaUI.stop_peak_training_timer!()
            @test HimalayaUI.PEAK_TRAINING.timer === nothing
            _pts_events!(db, 1)
            sleep(0.6)
            @test length(_pts_rows(db)) == 1
        end
    end
end

@testset "the boot starts the training timer only when it resumes work" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        try
            HimalayaUI.start_background!(db; resume = false)
            @test HimalayaUI.PEAK_TRAINING.timer === nothing
            HimalayaUI.stop_gc_timer!()
            HimalayaUI.start_background!(db; resume = false, train = true)
            @test HimalayaUI.PEAK_TRAINING.timer !== nothing
            HimalayaUI.stop_peak_training!()
            @test HimalayaUI.PEAK_TRAINING.timer === nothing
        finally
            HimalayaUI.stop_gc_timer!()
            HimalayaUI.stop_peak_training!()
            DBInterface.close!(db)
        end
    end
end

@testset "stop right after a start never waits for the run" begin
    mktempdir() do dir
        _pts_with_service(dir; command = _pts_fake("sleep(60)")) do db
            _pts_events!(db, 1)
            @test HimalayaUI.maybe_start_peak_training!(db; min_events = 1, quiet_seconds = 0.0) === :started
            t0 = time()
            HimalayaUI.stop_peak_training!()
            @test time() - t0 < 30
            @test HimalayaUI.PEAK_TRAINING.process === nothing && !HimalayaUI.peak_training_running()
            @test isempty(_pts_rows(db))
        end
    end
end

@testset "after a stop no check starts a child until the timer is started again" begin
    mktempdir() do dir
        _pts_with_service(dir; command = _pts_fake(_PTS_PRINT; env = ["FAKE_RESULT" => _pts_result(watermark = 1)])) do db
            _pts_events!(db, 1)
            kw = (min_events = 1, quiet_seconds = 0.0)
            HimalayaUI.stop_peak_training!()
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :stopped
            @test !HimalayaUI.peak_training_running() && HimalayaUI.PEAK_TRAINING.process === nothing
            @test isempty(_pts_rows(db))
            HimalayaUI.start_peak_training_timer!(db; interval = 3600, initial_delay = 3600, kw...)
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :started
            HimalayaUI.wait_peak_training()
            @test length(_pts_rows(db)) == 1
        end
    end
end

@testset "a run whose row cannot be written is not repeated until new labels arrive" begin
    mktempdir() do dir
        slow = "sleep(1.5); println(ENV[\"FAKE_RESULT\"])"
        _pts_with_service(dir; command = _pts_fake(slow; env = ["FAKE_RESULT" => _pts_result(watermark = 3)])) do db
            kw = (min_events = 3, quiet_seconds = 0.0)
            _pts_events!(db, 3)
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :started
            sleep(0.3)
            DBInterface.execute(db, "DROP TABLE peak_models")
            HimalayaUI.wait_peak_training()
            HimalayaUI.migrate_peak_models!(db)
            @test isempty(_pts_rows(db))
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :not_due
            _pts_events!(db, 3)
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :started
            HimalayaUI.wait_peak_training()
        end
    end
end

@testset "the default command caps the child's threads and lowers its priority" begin
    cmd = HimalayaUI.default_peak_training_command("/x/h.db"; extra = ["--trees", "5"])
    exec = cmd.exec
    threads = parse(Int, only(m.captures[1] for m in match.(r"^--threads=(\d+)$", exec) if m !== nothing))
    @test 1 <= threads <= max(1, Threads.nthreads(:default)) && threads <= max(1, Sys.CPU_THREADS ÷ 2)
    @test threads == HimalayaUI.peak_training_threads()
    if Sys.islinux() && Sys.which("nice") !== nothing
        @test basename(exec[1]) == "nice" && exec[2:3] == ["-n", "10"]
    end
    @test exec[end - 3:end] == ["--db", "/x/h.db", "--trees", "5"]
end

@testset "the real training command on a synthetic database records and promotes a run, and changes no reviewed peaks" begin
    mktempdir() do dir
        b = _pt_db(dir)
        db = b.db
        old = HimalayaUI.PEAK_TRAINING_COMMAND[]
        try
            _pts_events!(db, 1)
            snapshot() = (Tables.rowtable(DBInterface.execute(db, "SELECT exposure_id, q, found_by, model_id FROM auto_peaks ORDER BY id")),
                          Tables.rowtable(DBInterface.execute(db, "SELECT COUNT(*) AS n FROM user_actions")),
                          Tables.rowtable(DBInterface.execute(db, "SELECT id, analysis_version, reviewed_at FROM exposures ORDER BY id")))
            before = snapshot()
            HimalayaUI.PEAK_TRAINING.stopped = false
            HimalayaUI.PEAK_MODEL_PASS.stopped = false
            HimalayaUI.PEAK_TRAINING_COMMAND[] = path -> HimalayaUI.default_peak_training_command(path; extra = ["--trees", "20"])
            @test HimalayaUI.maybe_start_peak_training!(db; min_events = 1, quiet_seconds = 0.0) === :started
            HimalayaUI.wait_peak_training()
            HimalayaUI.wait_peak_model_pass()
            row = only(_pts_rows(db))
            @test row.status == "active"
            @test HimalayaUI.production_detector(db).id == row.id
            @test row.label_watermark == HimalayaUI.label_watermark(db)
            @test row.core_version == string(pkgversion(HimalayaUI.Himalaya))
            @test JSON3.read(row.metrics_json).gate.passed
            @test JSON3.read(row.counts_json).exposures == 36
            @test HPM.recipe_id(HPM.from_dict(HPM.Recipe, JSON3.read(row.recipe_json))) == row.recipe_id
            f = HPM.from_dict(HPM.Fitted, JSON3.read(row.fitted_json))
            @test HPM.canonical(HPM.to_dict(f)) == HPM.canonical(JSON3.read(row.fitted_json))
            @test length(readdir(joinpath(dir, "cache", "peak-features", "h"))) == 36
            @test isequal(snapshot(), before)
        finally
            HimalayaUI.PEAK_TRAINING_COMMAND[] = old
            HimalayaUI.stop_peak_training!()
            HimalayaUI.stop_peak_model_pass!()
            DBInterface.close!(db)
        end
    end
end

@testset "a stale model makes the service retrain once without new labels" begin
    mktempdir() do dir
        _pts_with_service(dir; command = _pts_fake(_PTS_PRINT; env = ["FAKE_RESULT" => _pts_result(status = "rejected", watermark = 1)])) do db
            HimalayaUI.PEAK_TRAINING.attempted_stale = false
            _pts_events!(db, 2)
            DBInterface.execute(db, """INSERT INTO peak_models (created_at, status, label_watermark, core_version, feature_version)
                                       VALUES ('t', 'retired', 2, '0', 'default-1;old')""")
            kw = (min_events = 20, quiet_seconds = 600.0)
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :started
            HimalayaUI.wait_peak_training()
            @test HimalayaUI.maybe_start_peak_training!(db; kw...) === :not_due
            @test length(_pts_rows(db)) == 2
        end
    end
end
