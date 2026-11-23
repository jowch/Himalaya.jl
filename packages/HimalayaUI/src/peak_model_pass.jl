mutable struct PeakModelPassState
    lock::ReentrantLock
    task::Union{Task, Nothing}
    running::Bool
    pending::Bool
    stopped::Bool
end

const PEAK_MODEL_PASS = PeakModelPassState(ReentrantLock(), nothing, false, false, false)

function peak_model_pass_ids(db::SQLite.DB, model_id::Integer)
    Int[Int(r.id) for r in sqlrows(DBInterface.execute(db, """
        SELECT id FROM exposures
         WHERE reviewed_at IS NULL AND trace_hash IS NOT NULL AND filename IS NOT NULL
           AND analysis_error IS NULL AND analysis_version = ?
           AND peak_model_version IS NOT ?
         ORDER BY experiment_id DESC, id""", [ANALYSIS_VERSION, Int(model_id)]))]
end

function run_peak_model_pass!(db::SQLite.DB)
    with_background_db(db) do bg
        pm = production_detector(bg)
        pm === nothing && return 0
        done = 0
        for id in peak_model_pass_ids(bg, pm.id)
            PEAK_MODEL_PASS.stopped && break
            current = production_detector(bg)
            (current === nothing || current.id != pm.id) && break
            try
                analyze_exposure!(bg, id; force_detect = true, only_unreviewed = true, log_event = false)
                done += 1
            catch err
                @warn "model pass could not re-detect exposure $id" exception = err
            end
            yield()
        end
        done
    end
end

function _peak_model_pass_loop(db::SQLite.DB)
    try
        while true
            try
                run_peak_model_pass!(db)
            catch err
                @warn "model pass failed" exception = err
            end
            again = lock(PEAK_MODEL_PASS.lock) do
                repeat = PEAK_MODEL_PASS.pending && !PEAK_MODEL_PASS.stopped
                PEAK_MODEL_PASS.pending = false
                PEAK_MODEL_PASS.running = repeat
                repeat
            end
            again || break
        end
    finally
        lock(() -> PEAK_MODEL_PASS.running = false, PEAK_MODEL_PASS.lock)
    end
    nothing
end

function start_peak_model_pass!(db::SQLite.DB)
    lock(PEAK_MODEL_PASS.lock) do
        PEAK_MODEL_PASS.stopped && return :stopped
        if PEAK_MODEL_PASS.running
            PEAK_MODEL_PASS.pending = true
            return :queued
        end
        PEAK_MODEL_PASS.running = true
        PEAK_MODEL_PASS.pending = false
        PEAK_MODEL_PASS.task = Threads.@spawn _peak_model_pass_loop(db)
        :started
    end
end

function wait_peak_model_pass()
    t = lock(() -> PEAK_MODEL_PASS.task, PEAK_MODEL_PASS.lock)
    t === nothing || wait(t)
    nothing
end

function stop_peak_model_pass!()
    t = lock(PEAK_MODEL_PASS.lock) do
        PEAK_MODEL_PASS.stopped = true
        PEAK_MODEL_PASS.task
    end
    t === nothing || wait(t)
    nothing
end

function boot_peak_model!(db::SQLite.DB)
    try
        retire_stale_peak_models!(db)
        active_peak_model_id(db) === nothing || start_peak_model_pass!(db)
    catch err
        @warn "peak model boot failed" exception = err
    end
    nothing
end
