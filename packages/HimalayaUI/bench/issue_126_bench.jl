# packages/HimalayaUI/bench/issue_126_bench.jl
#
# Microbench for issue #126. Run with:
#   julia --project=packages/HimalayaUI packages/HimalayaUI/bench/issue_126_bench.jl
#
# Measures `analyze_exposure!` cold and no-op latency on the example trace.
# Optional second arg adds N synthetic add-curations before the no-op runs
# to stress the curation-spam (replay-as-rerun) path.

using HimalayaUI
using HimalayaUI: open_db, create_schema!, create_experiment!, create_sample!,
                  create_exposure!, analyze_exposure!, load_dat
using SQLite, DBInterface
using Printf

const N_ADDS = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 0
const ITERS  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20

function setup_bench()
    tmp = mktempdir()
    analysis_dir = joinpath(tmp, "analysis")
    mkpath(analysis_dir)

    # Copy the test .dat into the analysis dir.
    src = joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat")
    dat_dest = joinpath(analysis_dir, "example_tot.dat")
    cp(src, dat_dest; force = true)

    db = open_db(joinpath(tmp, "himalaya.db"))
    create_schema!(db)
    exp_id = create_experiment!(db; path = tmp, data_dir = analysis_dir,
                                     analysis_dir = analysis_dir)
    s_id   = create_sample!(db; experiment_id = exp_id, name = "D1",
                                 display_name = "UX1")
    e_id   = create_exposure!(db; sample_id = s_id, filename = "example_tot")

    db, e_id, analysis_dir
end

function add_curations!(db, e_id, q_vals)
    for q in q_vals
        DBInterface.execute(db,
            "INSERT INTO peak_curations (exposure_id, q, kind) VALUES (?, ?, 'add')",
            [e_id, q])
    end
end

function bench(label, f, iters)
    f()  # warm up
    t0 = time_ns()
    for _ in 1:iters
        f()
    end
    elapsed_ms = (time_ns() - t0) / 1e6 / iters
    @printf "  %-30s %8.3f ms/iter\n" label elapsed_ms
    elapsed_ms
end

function main()
    db, e_id, analysis_dir = setup_bench()

    println("=== analyze_exposure! microbench ===")
    println("iters=$ITERS  n_add_curations=$N_ADDS")

    # Cold run (no prior state).
    print("Cold: ")
    t0 = time_ns()
    analyze_exposure!(db, e_id, analysis_dir)
    @printf "%.3f ms\n" ((time_ns() - t0) / 1e6)

    if N_ADDS > 0
        # Synthetic add curations near observed peak q's.
        q_arr, _, _ = load_dat(joinpath(analysis_dir, "example_tot.dat"))
        step = max(1, length(q_arr) ÷ (N_ADDS + 1))
        add_qs = [q_arr[i * step] for i in 1:N_ADDS]
        add_curations!(db, e_id, add_qs)
        # Re-run once so the inputs_hash and indices reflect the curation set.
        analyze_exposure!(db, e_id, analysis_dir)
    end

    # No-op path: trace unchanged, indices already exist.
    bench("no-op (trace_known)", () -> analyze_exposure!(db, e_id, analysis_dir;
                                                         trace_known_unchanged = true),
          ITERS)
    bench("no-op (trace_unknown)", () -> analyze_exposure!(db, e_id, analysis_dir),
          ITERS)
end

main()
