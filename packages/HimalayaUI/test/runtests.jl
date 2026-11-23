using Test

const GROUP = get(ENV, "GROUP", "All")
const _TEST_TIMES = Dict{String,Float64}()
_timed_include(f) = (_TEST_TIMES[f] = @elapsed include(f))
_want(g) = GROUP == "All" || GROUP == g

# Buckets balanced from M0.2 timing (heaviest files spread across buckets).
# IMPORTANT: under GROUP=All these run in the SAME total order as before
# (db → pipeline → routes → events → wire → …). Each bucket includes
# test_http.jl first (the isdefined guard makes re-include a no-op).
# Every bucket lists test_http.jl FIRST — it transitively includes test_fixtures.jl
# (M3.0), test_inproc.jl, test_template_db.jl, so any cross-file helper is in scope
# regardless of which bucket a file lands in. (test_http.jl defines no tests, so the
# Pass-sum invariant across buckets holds. The seen-set dedup keeps it included once
# per process.)
const GROUPS = [
    ("db",       ["test_http.jl","test_config.jl","test_db.jl",
                  "test_migrate_speculative_durability.jl","test_migrate_hex_sqrt11.jl",
                  "test_migrate_ratio_position_int.jl",
                  "test_ingestion_schema.jl"]),
    ("pipeline", ["test_http.jl","test_datfile.jl","test_hash.jl","test_hash_peak_set_memoization.jl",
                  "test_pipeline.jl","test_auto_group_peak_id_claiming.jl",
                  "test_effective_peaks_sharpness_passthrough.jl","test_ingestion_core.jl",
                  "test_ingest_reconcile.jl","test_setups.jl","test_beamtime.jl","test_routes_setups.jl","test_json.jl","test_image.jl","test_demo_data.jl",
                  "test_peak_labels.jl"]),
    ("routes",   ["test_http.jl","test_routes_experiments.jl","test_experiments_patch.jl",
                  "test_routes_samples.jl","test_routes_exposures.jl","test_routes_image.jl",
                  "test_routes_peaks.jl","test_peak_matching.jl","test_routes_trace.jl","test_routes_export.jl",
                  "test_route_response_shapes.jl",
                  "test_routes_fs.jl"]),
    ("routes2",  ["test_http.jl","test_routes_users.jl","test_ingestion_scan_api.jl",
                  "test_routes_analysis.jl","test_speculative.jl",
                  "test_route_validation_routing.jl",
                  "test_inproc_equivalence.jl","test_routes_bridges.jl","test_bridge_calls.jl","test_routes_compares.jl","test_routes_labels.jl","test_routes_label_colors.jl","test_label_suggestions.jl"]),
    ("events",   ["test_http.jl","test_actions.jl","test_events.jl","test_assignment_reattach.jl",
                  "test_ingestion_structural_events.jl",
                  "test_assignments.jl","test_sample_calls.jl","test_fast_skip.jl","test_idempotency.jl",
                  "test_idempotency_replay_invariant.jl","test_concurrent_writes.jl","test_db_connections.jl",
                  "test_idempotency_sse_suppression.jl","test_peak_provenance.jl"]),
    ("wire",     ["test_http.jl","test_routes_health.jl","test_routes_status.jl",
                  "test_sse.jl","test_routes_sse_broadcast.jl","test_spa_fallback.jl",
                  "test_precompile_workload.jl"]),
    ("migration", ["test_http.jl","test_migration_upgrade.jl"]),
]

GROUP == "All" || GROUP in first.(GROUPS) ||
    error("runtests: unknown GROUP=$(GROUP); expected All or one of $(join(first.(GROUPS), ", "))")

# Single source of truth for the serial bisect order — the exact historical
# runtests include order. The buckets above must partition this same set.
const ALL_ORDER = ["test_config.jl","test_db.jl",
                   "test_migrate_speculative_durability.jl","test_migrate_hex_sqrt11.jl",
                  "test_migrate_ratio_position_int.jl",
                   "test_ingestion_schema.jl",
                   "test_datfile.jl","test_hash.jl",
                   "test_hash_peak_set_memoization.jl","test_pipeline.jl",
                   "test_auto_group_peak_id_claiming.jl",
                   "test_effective_peaks_sharpness_passthrough.jl","test_ingestion_core.jl",
                   "test_ingest_reconcile.jl","test_setups.jl","test_beamtime.jl","test_routes_setups.jl","test_fast_skip.jl",
                   "test_json.jl","test_demo_data.jl","test_peak_labels.jl","test_http.jl","test_inproc_equivalence.jl",
                   "test_routes_health.jl","test_routes_users.jl","test_routes_experiments.jl",
                   "test_experiments_patch.jl",
                   "test_ingestion_scan_api.jl",
                   "test_ingestion_structural_events.jl",
                   "test_routes_samples.jl","test_routes_exposures.jl","test_image.jl",
                   "test_routes_image.jl","test_routes_status.jl","test_routes_peaks.jl","test_peak_matching.jl",
                   "test_routes_trace.jl","test_routes_analysis.jl",
                   "test_speculative.jl","test_routes_export.jl",
                   "test_actions.jl","test_events.jl","test_assignments.jl","test_sample_calls.jl",
                   "test_assignment_reattach.jl","test_sse.jl","test_routes_sse_broadcast.jl",
                   "test_route_response_shapes.jl","test_route_validation_routing.jl",
                   "test_idempotency.jl","test_idempotency_replay_invariant.jl",
                   "test_concurrent_writes.jl","test_db_connections.jl","test_idempotency_sse_suppression.jl",
                   "test_routes_bridges.jl","test_bridge_calls.jl","test_routes_compares.jl","test_routes_labels.jl","test_routes_label_colors.jl","test_label_suggestions.jl",
                   "test_routes_fs.jl",
                   "test_peak_provenance.jl","test_spa_fallback.jl","test_precompile_workload.jl","test_migration_upgrade.jl"]

# Drift guard: the buckets and ALL_ORDER must cover the identical file set, so a
# parallel / GROUP=<name> run can never silently skip (or double-run) a file
# relative to GROUP=All. Fails loudly at load if a test file is added to one list
# but not the other. (test_http.jl is in every bucket + ALL_ORDER → in both sets.)
let bucket_files = Set(Iterators.flatten(fs for (_, fs) in GROUPS))
    drift = symdiff(bucket_files, Set(ALL_ORDER))
    isempty(drift) || error("runtests: GROUP buckets ↔ ALL_ORDER drift: $(sort(collect(drift)))")
    listed = [f for (_, fs) in GROUPS for f in fs if f != "test_http.jl"]
    dupes = sort(unique(f for f in listed if count(==(f), listed) > 1))
    isempty(dupes) || error("runtests: test files in more than one GROUP bucket: $(dupes)")
end

# Second guard (closes the disk↔list gap): every test_*.jl on disk must be in
# ALL_ORDER, so a newly-added file can't silently run in NO configuration.
# Exclude the include-only helpers — they define no testsets and are pulled in
# transitively via test_http.jl, so they are intentionally absent from ALL_ORDER.
let helpers  = Set(["test_http.jl","test_fixtures.jl","test_inproc.jl","test_template_db.jl"]),
    on_disk  = Set(f for f in readdir(@__DIR__)
                   if startswith(f, "test_") && endswith(f, ".jl") && !(f in helpers)),
    orphaned = setdiff(on_disk, Set(ALL_ORDER))
    isempty(orphaned) || error("runtests: test files on disk but in no run list: $(sort(collect(orphaned)))")
end

# GROUP=All must reproduce the historical order exactly — assert coverage.
@testset "HimalayaUI" begin
    if GROUP == "All"
        for f in ALL_ORDER
            _timed_include(f)
        end
    else
        seen = Set{String}()
        for (name, files) in GROUPS
            _want(name) || continue
            for f in files
                f in seen && continue   # test_http.jl appears in several buckets
                push!(seen, f)
                _timed_include(f)
            end
        end
    end
end

let rows = sort(collect(_TEST_TIMES); by = last, rev = true)
    println("\n── per-file test timing (slowest first) ──")
    for (f, t) in rows
        println(rpad(f, 42), lpad(round(t; digits=1), 8), " s")
    end
end
