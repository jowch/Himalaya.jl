using HimalayaUI, HTTP, JSON3, Oxygen, SQLite

const H = HimalayaUI

const BUDGET_MS = Dict("boot" => 200, "experiment" => 200, "first sample" => 150)
budget(step) = get(BUDGET_MS, step, 100) * parse(Float64, get(ENV, "FIRST_USE_BUDGET_SCALE", "1"))

const SERVER_COMPILE_NS = Ref(UInt64(0))
const TRACE = Ref{Union{Nothing, String}}(nothing)
const SERVER_TRACE = String[]

trace_length() = (Base.Libc.flush_cstdio(); countlines(TRACE[]))

function on_server(f)
    from = TRACE[] === nothing ? 0 : trace_length()
    t0 = Base.cumulative_compile_time_ns()[1]
    out = f()
    SERVER_COMPILE_NS[] += Base.cumulative_compile_time_ns()[1] - t0
    TRACE[] === nothing || append!(SERVER_TRACE, readlines(TRACE[])[from+1:trace_length()])
    out
end

function call(method, target, body = nothing)
    payload = body === nothing ? UInt8[] : Vector{UInt8}(JSON3.write(body))
    res = on_server(() -> H._workload_inproc(method, target, payload))
    res.status < 300 || error("$method $target returned $(res.status): $(String(copy(res.body)))")
    isempty(res.body) || !occursin("json", HTTP.header(res, "Content-Type")) ? nothing : JSON3.read(res.body)
end

function prepare(root, db_path)
    db = H.open_db(db_path)
    H.boot!(db)
    call("POST", "/api/users", Dict(:username => "precompile"))
    rs = only(call("GET", "/api/fs/preview?path=" * HTTP.escapeuri(root)).setups)
    call("POST", "/api/experiments", Dict(:name => "Demo (synthetic)", :path => root,
        :setups => [Dict(:root => rs.root, :data_dir => rs.data_dir, :analysis_dir => rs.analysis_dir,
                         :included => true, :patterns => rs.patterns, :geometry => rs.geometry)]))
    deadline = time() + 600
    while time() < deadline
        e = only(call("GET", "/api/experiments"))
        e = call("GET", "/api/experiments/$(e.id)")
        e.ingest_status == "complete" && return println("Ingested experiment $(e.id).")
        e.ingest_status == "failed" && error("ingest failed: $(get(e, :ingest_error, ""))")
        sleep(0.5)
    end
    error("ingest did not finish within 10 minutes")
end

function sample_open(sample)
    exposures = call("GET", "/api/samples/$(sample.id)/exposures")
    x = sample.analyzed_exposure_id
    call("GET", "/api/samples/$(sample.id)/call")
    for e in exposures
        get(e, :kind, "file") == "file" && call("GET", "/api/exposures/$(e.id)/image?thumb=1")
    end
    for part in ("trace", "peaks", "indices", "assignment")
        call("GET", "/api/exposures/$x/$part")
    end
    images = [e for e in exposures if get(e, :kind, "file") == "file" && get(e, :selected, false) === true]
    call("GET", "/api/exposures/$(first(images).id)/image")
    length(images) > 1 && foreach(e -> call("GET", "/api/exposures/$(e.id)/trace"), images)
    x
end

function measure(db_path)
    trace = Base.JLOptions().trace_compile
    trace == C_NULL && error("run with --trace-compile=<file> --trace-compile-timing, so a failure can name what compiled")
    db = H.open_db(db_path)
    rm(dirname(H.thumb_cache_dir(db)); recursive = true, force = true)
    TRACE[] = unsafe_string(trace)
    on_server(() -> nothing)
    SERVER_COMPILE_NS[] = 0
    failures = String[]
    results = Pair{String, Float64}[]
    function step(f, name)
        empty!(SERVER_TRACE)
        t0 = SERVER_COMPILE_NS[]
        out = f()
        ms = (SERVER_COMPILE_NS[] - t0) / 1e6
        push!(results, name => ms)
        if ms > budget(name)
            timed = sort!([(parse(Float64, m.captures[1]), l) for l in SERVER_TRACE
                           for m in (match(r"^#=\s*([\d.]+) ms", l),) if m !== nothing]; rev = true)
            push!(failures, "$name compiled for $(round(Int, ms)) ms (budget $(round(Int, budget(name))) ms). Slowest:\n" *
                  join(("    " * first(l, 300) for (_, l) in first(timed, 8)), "\n"))
        end
        out
    end

    step("boot") do
        on_server() do
            wait(H.boot!(db; resume = false))
            H.stop_gc_timer!()
        end
    end
    experiment = step("experiment") do
        call("GET", "/api/users"); exps = call("GET", "/api/experiments"); e = only(exps).id
        call("GET", "/api/experiments/$e")
        for part in ("bridges", "problems", "compares", "labels")
            call("GET", "/api/experiments/$e/$part")
        end
        call("GET", "/api/label-colors")
        call("GET", "/api/samples?experiment_id=$e")
        e
    end
    samples = [s for s in call("GET", "/api/samples?experiment_id=$experiment")
               if s.analyzed_exposure_id !== nothing]
    length(samples) >= 3 || error("the demo needs at least 3 samples with a kept image")
    step(() -> sample_open(samples[1]), "first sample")
    step(() -> sample_open(samples[2]), "second sample")

    target = findfirst(smp -> any(p -> p.source == "auto",
                                  call("GET", "/api/exposures/$(smp.analyzed_exposure_id)/peaks")), samples)
    target === nothing && error("no demo sample has auto peaks on its shown image")
    x = samples[target].analyzed_exposure_id
    s = samples[target].id
    added = step(() -> call("POST", "/api/exposures/$x/peaks", Dict(:q => 0.123)), "add peak")
    step(() -> call("DELETE", "/api/peaks/manual/$(added.id)"), "remove peak")
    auto = first(p for p in call("GET", "/api/exposures/$x/peaks") if p.source == "auto")
    step(() -> call("PATCH", "/api/peaks/auto/$(auto.id)", Dict(:excluded => true)), "exclude peak")
    step(() -> call("PATCH", "/api/peaks/auto/$(auto.id)", Dict(:excluded => false)), "restore peak")
    indices = call("GET", "/api/exposures/$x/indices")
    if !isempty(indices)
        i = first(indices).id
        step(() -> call("POST", "/api/exposures/$x/assignment/members", Dict(:index_id => i)), "tick index")
        step(() -> call("DELETE", "/api/exposures/$x/assignment/members/$i"), "untick index")
    end
    step(() -> call("POST", "/api/samples/$s/labels", Dict(:label => "Buffer")), "add label")
    step(() -> call("DELETE", "/api/samples/$s/labels?label=Buffer"), "remove label")
    other = samples[mod1(target + 1, length(samples))].id
    step(() -> call("PATCH", "/api/samples/$s/compares", Dict(:add => [other])), "add comparison")
    step(() -> call("PATCH", "/api/samples/$s/compares", Dict(:remove => [other])), "remove comparison")
    step(() -> call("POST", "/api/exposures/$x/custom-index", Dict(:phase => "Ia3d", :basis => 0.05)), "custom index")
    images = [e.id for e in call("GET", "/api/samples/$s/exposures") if e.kind == "file" && e.selected === true]
    step(() -> call("PUT", "/api/samples/$s/selection", Dict(:exposure_ids => images[1:1])), "select one image")
    step(() -> call("PUT", "/api/samples/$s/selection", Dict(:exposure_ids => images)), "select images")
    img = first(images)
    step(() -> call("PATCH", "/api/exposures/$img/status", Dict(:status => "rejected")), "drop image")
    step(() -> call("PATCH", "/api/exposures/$img/status", Dict(:status => nothing)), "keep image")

    for (name, ms) in results
        println(rpad(name, 20), lpad(round(Int, ms), 6), " ms compiling   (budget ", round(Int, budget(name)), " ms)")
    end
    if !isempty(failures)
        println(stderr, """

            First use compiled more than its budget. Each step below is the first time a fresh server \
            does it on the demo. Usual causes: a route the precompile workload does not drive, a value \
            whose type depends on the data (read DB rows with `sqlrows`), or cached code thrown away at \
            load (run scripts/check_precompile_cache.jl). See packages/HimalayaUI/src/AGENTS.md, "Cold start".
            """)
        foreach(f -> println(stderr, f, "\n"), failures)
        exit(1)
    end
end

if length(ARGS) == 3 && ARGS[1] == "prepare"
    prepare(abspath(ARGS[2]), ARGS[3])
elseif length(ARGS) == 2 && ARGS[1] == "measure"
    Base.cumulative_compile_timing(true)
    measure(ARGS[2])
else
    println(stderr, "usage: check_first_use.jl prepare <demo root> <db>  |  check_first_use.jl measure <db>")
    exit(2)
end
