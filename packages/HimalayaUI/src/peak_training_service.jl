using JSON3, SQLite, DBInterface

const PEAK_TRAINING_CHECK_SECONDS = 300.0
const PEAK_TRAINING_FIRST_CHECK_SECONDS = 60.0
const PEAK_TRAINING_STDERR_TAIL = 1700

mutable struct PeakTrainingState
    lock::ReentrantLock
    task::Union{Task, Nothing}
    process::Union{Base.Process, Nothing}
    cancelled::Bool
    timer::Union{Timer, Nothing}
    exit_hook::Bool
    stopped::Bool
    attempted_watermark::Int
end

const PEAK_TRAINING = PeakTrainingState(ReentrantLock(), nothing, nothing, false, nothing, false, false, 0)

peak_training_threads() = max(1, min(Threads.nthreads(:default), Sys.CPU_THREADS ÷ 2))

function default_peak_training_command(db_path::AbstractString; extra::Vector{String} = String[])
    cmd = `$(Base.julia_cmd()) --threads=$(peak_training_threads()) --project=$(dirname(Base.active_project())) -e "using HimalayaUI; main(ARGS)" -- train-peak-model --db $db_path $extra`
    nice = Sys.islinux() ? Sys.which("nice") : nothing
    nice === nothing ? cmd : `$nice -n 10 $cmd`
end

const PEAK_TRAINING_COMMAND = Ref{Function}(default_peak_training_command)

peak_training_running() = lock(PEAK_TRAINING.lock) do
    t = PEAK_TRAINING.task
    t !== nothing && !istaskdone(t)
end

function _stderr_tail(path)
    isfile(path) || return ""
    s = read(path, String)
    lastindex(s) <= PEAK_TRAINING_STDERR_TAIL ? s : s[prevind(s, lastindex(s) - PEAK_TRAINING_STDERR_TAIL + 1):end]
end

function _record_training_outcome!(db::SQLite.DB, watermark::Integer, outcome)
    with_background_db(db) do bg
        lock(_DB_WRITE_LOCK) do
            if outcome isa AbstractString
                insert_failed_peak_model!(bg; watermark, error = outcome)
            else
                insert_peak_model!(bg, outcome)
            end
        end
    end
end

function _run_peak_training_child(db::SQLite.DB, watermark::Integer)
    dir = mktempdir()
    try
        out = joinpath(dir, "result.json"); err = joinpath(dir, "stderr.txt")
        cmd = PEAK_TRAINING_COMMAND[](String(db.file))
        proc = lock(PEAK_TRAINING.lock) do
            PEAK_TRAINING.cancelled && return nothing
            p = run(pipeline(ignorestatus(cmd); stdout = out, stderr = err); wait = false)
            PEAK_TRAINING.process = p
            p
        end
        proc === nothing && return nothing
        wait(proc)
        PEAK_TRAINING.cancelled && return nothing
        outcome = if !success(proc)
            "peak model training exited with code $(proc.exitcode): $(_stderr_tail(err))"
        else
            try
                r = JSON3.read(read(out, String))
                Dict{String, Any}(
                    "status" => String(r.status), "watermark" => Int(r.watermark),
                    "counts" => r.counts, "recipe" => r.recipe, "recipe_id" => r.recipe_id,
                    "feature_version" => r.feature_version, "core_version" => String(r.core_version),
                    "threshold" => r.threshold, "nms" => r.nms, "metrics" => r.metrics, "fitted" => r.fitted)
            catch e
                "peak model training printed no usable result ($(first(sprint(showerror, e), 150))): $(_stderr_tail(err))"
            end
        end
        _record_training_outcome!(db, outcome isa AbstractString ? watermark : outcome["watermark"], outcome)
        if outcome isa AbstractString
            @warn "peak model training failed" error = outcome
        else
            @info "peak model training finished" status = outcome["status"] watermark = outcome["watermark"] threshold = outcome["threshold"] reasons = collect(outcome["metrics"].gate.reasons)
        end
        nothing
    catch e
        PEAK_TRAINING.cancelled || @warn "peak model training run failed" exception = (e, catch_backtrace())
        nothing
    finally
        lock(() -> PEAK_TRAINING.process = nothing, PEAK_TRAINING.lock)
        rm(dir; recursive = true, force = true)
    end
end

function maybe_start_peak_training!(db::SQLite.DB; min_events::Integer = PEAK_TRAINING_MIN_EVENTS,
                                    quiet_seconds::Real = PEAK_TRAINING_QUIET_SECONDS)
    _poolable(db.file) || return :unavailable
    lock(PEAK_TRAINING.lock) do
        PEAK_TRAINING.stopped && return :stopped
        due = with_background_db(bg -> peak_training_due(bg; min_events, quiet_seconds,
                                                         after = PEAK_TRAINING.attempted_watermark), db)
        due || return :not_due
        peak_training_running() && return :busy
        watermark = with_background_db(label_watermark, db)
        PEAK_TRAINING.attempted_watermark = max(PEAK_TRAINING.attempted_watermark, watermark)
        PEAK_TRAINING.cancelled = false
        PEAK_TRAINING.task = Threads.@spawn _run_peak_training_child(db, watermark)
        :started
    end
end

function wait_peak_training()
    t = lock(() -> PEAK_TRAINING.task, PEAK_TRAINING.lock)
    t === nothing || wait(t)
    nothing
end

function stop_peak_training!()
    stop_peak_training_timer!()
    task, proc = lock(PEAK_TRAINING.lock) do
        PEAK_TRAINING.stopped = true
        PEAK_TRAINING.cancelled = true
        PEAK_TRAINING.task, PEAK_TRAINING.process
    end
    proc === nothing || process_running(proc) && kill(proc)
    task === nothing || wait(task)
    nothing
end

function start_peak_training_timer!(db::SQLite.DB; interval::Real = PEAK_TRAINING_CHECK_SECONDS,
                                    initial_delay::Real = PEAK_TRAINING_FIRST_CHECK_SECONDS, kwargs...)
    stop_peak_training_timer!()
    lock(PEAK_TRAINING.lock) do
        PEAK_TRAINING.stopped = false
        PEAK_TRAINING.attempted_watermark = 0
    end
    if !PEAK_TRAINING.exit_hook
        PEAK_TRAINING.exit_hook = true
        atexit(stop_peak_training!)
    end
    PEAK_TRAINING.timer = Timer(initial_delay; interval) do _
        try
            maybe_start_peak_training!(db; kwargs...)
        catch err
            @warn "peak model training check failed" exception = err
        end
    end
    nothing
end

function stop_peak_training_timer!()
    t = PEAK_TRAINING.timer
    PEAK_TRAINING.timer = nothing
    t === nothing || close(t)
    nothing
end
