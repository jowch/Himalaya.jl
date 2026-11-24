# Ingest the synthetic demo experiment into a RUNNING HimalayaUI backend, through
# the same two calls the New-experiment screen makes:
#   1. GET  /api/fs/preview   (layout, counts and geometry for each setup folder)
#   2. POST /api/experiments  (create + async first scan)
# then waits for the scan to finish. `make demo` runs this; `make demo-empty` skips
# it so the funnel itself can be walked in the browser.
#
# Usage:
#   julia --project=packages/HimalayaUI packages/HimalayaUI/scripts/demo_ingest.jl \
#       <experiment_root>... [--port 8080]
# Several roots are ingested one after another (`make demo-large`); each experiment
# is named after its directory then.

using HTTP, JSON3

function main(args = ARGS)
    roots = String[]; port = 8080
    i = 1
    while i <= length(args)
        if args[i] == "--port"
            i < length(args) || error("--port needs a value")
            port = parse(Int, args[i+1]); i += 1
        elseif startswith(args[i], "-")
            error("unknown option $(args[i]) (usage: demo_ingest.jl [root] [--port N])")
        else
            push!(roots, args[i])
        end
        i += 1
    end
    isempty(roots) && push!(roots, "scratch/demo-experiment")
    names = length(roots) == 1 ? ["Demo (synthetic)"] :
        ["Demo " * basename(rstrip(abspath(r), '/')) for r in roots]
    allunique(names) || error("experiment names collide: $names")
    for (root, name) in zip(roots, names)
        ingest_root(abspath(root), port, name)
    end
end

function ingest_root(root, port, name)
    base = "http://127.0.0.1:$port/api"
    getj(path, q) = JSON3.read(HTTP.get(base * path; query = q).body)

    setups = sort(collect(getj("/fs/preview", ["path" => root]).setups); by = rs -> rs.mode == "autosampler")
    isempty(setups) && error("$root holds no setup folder. Regenerate it with `make demo-clean demo-data`.")
    if !all(rs -> rs.layout_detected, setups)
        error("$root doesn't look like an SSRL-style demo tree (no analysis/.../tot_files " *
              "*_tot.dat found). Regenerate it with `make demo-clean demo-data`.")
    end
    data_dirs = Set(rs.data_dir for rs in setups)

    # Already there? Done if its scan finished; otherwise (a previous run was
    # interrupted mid-scan, or the scan failed) force a rescan and wait for it.
    for e in JSON3.read(HTTP.get("$base/experiments").body)
        any(st -> get(st, :data_dir, nothing) in data_dirs, get(e, :setups, [])) || continue
        if get(e, :ingest_status, "") == "complete"
            println("$name already ingested (id $(e.id)); nothing to do.")
            return
        end
        println("Demo experiment $(e.id) is '$(get(e, :ingest_status, ""))'; rescanning.")
        HTTP.post("$base/experiments/$(e.id)/scan", ["Content-Type" => "application/json"],
                  JSON3.write(Dict(:force => true)))
        return wait_for_scan(base, e.id)
    end

    println("Found $(length(setups)) setups, $(sum(rs -> rs.exposures, setups)) exposures, ",
            "$(sum(rs -> rs.without_trace, setups)) without a trace.")

    body = Dict(
        :name => name, :path => root,
        :setups => [Dict(:root => rs.root, :data_dir => rs.data_dir, :analysis_dir => rs.analysis_dir,
                         :included => true, :patterns => rs.patterns, :geometry => rs.geometry)
                    for rs in setups],
    )
    created = JSON3.read(HTTP.post("$base/experiments", ["Content-Type" => "application/json"],
                                   JSON3.write(body)).body)
    println("Created experiment $(created.id).")
    wait_for_scan(base, created.id)
end

function wait_for_scan(base, id)
    print("Scanning")
    for _ in 1:600
        st = JSON3.read(HTTP.get("$base/experiments/$id").body)
        status = get(st, :ingest_status, "")
        if status == "complete"
            println(" done. ", get(st, :stats, ""))
            return
        elseif status == "failed"
            println(); error("scan failed; see the backend log")
        end
        print("."); sleep(1)
    end
    println(); error("timed out waiting for the scan")
end

main()
