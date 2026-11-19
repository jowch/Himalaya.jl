# Ingest the synthetic demo experiment into a RUNNING HimalayaUI backend, through
# the same three calls the New-experiment funnel makes:
#   1. GET  /api/fs/resolve   (layout auto-detect: data_dir, analysis_dir, patterns, setup file)
#   2. GET  /api/fs/manifest  (file match + geometry derivation with that setup file)
#   3. POST /api/experiments  (create + async first scan)
# then waits for the scan to finish. `make demo` runs this; `make demo-empty` skips
# it so the funnel itself can be walked in the browser.
#
# Usage:
#   julia --project=packages/HimalayaUI packages/HimalayaUI/scripts/demo_ingest.jl \
#       <experiment_root> [--port 8080]

using HTTP, JSON3

function main(args = ARGS)
    root = abspath(get(args, 1, "scratch/demo-experiment"))
    port = (i = findfirst(==("--port"), args)) === nothing ? 8080 : parse(Int, args[i+1])
    base = "http://127.0.0.1:$port/api"
    getj(path, q) = JSON3.read(HTTP.get(base * path; query = q).body)

    existing = JSON3.read(HTTP.get("$base/experiments").body)
    rs = getj("/fs/resolve", ["path" => root])
    for e in existing
        if get(e, :data_dir, nothing) == rs.data_dir
            println("Demo experiment already ingested (id $(e.id)); nothing to do.")
            return
        end
    end

    mf = getj("/fs/manifest", [
        "path" => rs.data_dir, "analysis_dir" => rs.analysis_dir,
        "image_pattern" => rs.image_pattern, "metadata_pattern" => rs.metadata_pattern,
        "integration_pattern" => rs.integration_pattern,
        "setup_file" => something(rs.setup_file, ""),
    ])
    println("Matched: ", mf.matched)

    body = Dict(
        :name => "Demo (synthetic)", :data_dir => rs.data_dir, :analysis_dir => rs.analysis_dir,
        :patterns => Dict(:image => rs.image_pattern, :metadata => rs.metadata_pattern,
                          :integration => rs.integration_pattern),
        :geometry => get(mf, :geometry, Dict()),
    )
    created = JSON3.read(HTTP.post("$base/experiments", ["Content-Type" => "application/json"],
                                   JSON3.write(body)).body)
    print("Created experiment $(created.id); scanning")
    for _ in 1:600
        st = JSON3.read(HTTP.get("$base/experiments/$(created.id)").body)
        status = get(st, :ingest_status, "")
        if status == "complete"
            println(" done. ", get(st, :stats, ""))
            return
        elseif status == "failed"
            println(); error("first scan failed; see the backend log")
        end
        print("."); sleep(1)
    end
    println(); error("timed out waiting for the first scan")
end

main()
