# packages/HimalayaUI/src/grouping.jl
#
# Directory scanning + sample grouping (spec §5, §9.1).
#
# scan_directory: enumerate (tif, prp, dat) triplets from a beamtime directory.
# group_into_samples: cluster exposures into loads and samples (Tasks 5–7).
# scan_and_group!: transactional orchestrator (ingest.jl, Task 8).
#
# Naming note: `auto_group` already exists in pipeline.jl for index-candidate
# grouping (an unrelated concept). This module uses `group_into_samples` (spec §9.1).

"""
    ExposureMeta

Per-file metadata collected during directory scan. All fields except `stem`
are nullable (`nothing` or `missing`).
"""
struct ExposureMeta
    stem         ::String
    tif_path     ::Union{String, Nothing}
    prp_path     ::Union{String, Nothing}
    dat_path     ::Union{String, Nothing}
    prp          ::Union{NamedTuple, Nothing}  # parse_prp result; nothing if no .prp
end

Base.haskey(m::ExposureMeta, k::Symbol) = k in fieldnames(ExposureMeta)

"""
    scan_directory(data_dir, analysis_dir;
                   tif_pattern = "{name}.tif",
                   prp_pattern = "{name}.prp",
                   dat_pattern = "{name}.dat") -> Vector{ExposureMeta}

Enumerate every TIFF found in `data_dir`, then pair each with its PRP and .dat
sidecar. Returns one `ExposureMeta` per TIFF stem, sorted by stem.

Reuses the `resolve_file_path` logic from `config.jl` for sidecar lookup once
a stem is known, constructing a minimal `ExperimentConfig` for dispatch (resolves
the §9.1 open question: loosen dispatch vs. construct minimal config — we construct
a minimal config since `resolve_file_path` takes `ExperimentConfig` for dispatch
only and none of its fields are read by the dispatch branch).

TIF enumeration does NOT use `resolve_files` with an empty prefix because
`_matches_prefix_with_boundary` rejects filenames whose first character is
alphanumeric when the prefix is empty — which is the common case. Instead we
scan the directory directly and strip the TIF suffix to derive stems.
"""
function scan_directory(
    data_dir::AbstractString,
    analysis_dir::AbstractString;
    tif_pattern::String = "{name}.tif",
    prp_pattern::String = "{name}.prp",
    dat_pattern::String = "{name}.dat",
)::Vector{ExposureMeta}
    # Build a minimal ExperimentConfig for the dispatch arg to resolve_file_path.
    # The dispatch method only needs the struct type; no fields are read.
    cfg = _minimal_scan_config(tif_pattern, prp_pattern, dat_pattern,
                               String(data_dir), String(analysis_dir))

    # Derive the literal TIF suffix from the pattern (e.g. "{name}.tif" → ".tif").
    tif_parts = split(tif_pattern, "{name}"; limit=2)
    length(tif_parts) == 2 || error("tif_pattern must contain {{name}}: $tif_pattern")
    tif_prefix_literal = String(tif_parts[1])  # e.g. "" for "{name}.tif"
    tif_suffix = String(tif_parts[2])           # e.g. ".tif"

    # Enumerate all files in data_dir that match the TIF pattern (exact suffix).
    scan_subdir = dirname(tif_prefix_literal)
    tif_scan_dir = isempty(scan_subdir) ? String(data_dir) : joinpath(data_dir, scan_subdir)

    if !isdir(tif_scan_dir)
        return ExposureMeta[]
    end

    tif_files = filter(readdir(tif_scan_dir)) do f
        startswith(f, tif_prefix_literal) && endswith(f, tif_suffix)
    end
    sort!(tif_files)

    # Strip the TIF suffix to get stems, then pair each stem with its sidecars.
    suffix_len = length(tif_suffix)
    metas = ExposureMeta[]
    for fname in tif_files
        stem = fname[1:end - suffix_len]
        tif_path  = joinpath(tif_scan_dir, fname)
        prp_path  = resolve_file_path(cfg, data_dir, stem, prp_pattern)
        dat_path  = resolve_file_path(cfg, analysis_dir, stem, dat_pattern)
        prp_parsed = prp_path !== nothing ? parse_prp(prp_path) : nothing
        push!(metas, ExposureMeta(stem, tif_path, prp_path, dat_path, prp_parsed))
    end
    return metas
end

"""
    _minimal_scan_config(tif_pattern, prp_pattern, dat_pattern, data_dir, analysis_dir) -> ExperimentConfig

Construct the minimal `ExperimentConfig` needed to dispatch `resolve_files` /
`resolve_file_path`. The manifest and beamline fields are set to safe sentinel
values and are never read during a scan.

Field order confirmed from config.jl struct definition (23 fields total):
  name, description, manifest_file                         # [experiment]
  energy_kev, flight_path_m, q_units,
    beam_center_x, beam_center_y, pixel_size_um           # [beamline]
  delimiter, skip_rows, header_row,
    col_sample_id, col_name, col_display_name,
    col_filenames, col_notes_sample, col_notes_exposure    # [manifest]
  data_dir, analysis_dir, exposure_type                    # [layout]
  integration_pattern, image_pattern                       # [files]

If ExperimentConfig gains new fields, this positional call must be updated
in lockstep — re-read config.jl before editing.
"""
# ---------------------------------------------------------------------------
# Time-gap load segmentation (spec §5, step 2)
# ---------------------------------------------------------------------------

"""
    _segment_loads(metas; gap_k = 10.0) -> Vector{Vector{ExposureMeta}}

Split `metas` (sorted by timestamp) into loads using gap-relative segmentation:
a time gap > `gap_k × median(all_consecutive_gaps)` starts a new load.

Returns only the load groups (no flag). Use `_segment_loads_with_flag` to
distinguish the unimodal fallback.
"""
function _segment_loads(metas::Vector{ExposureMeta}; gap_k::Float64 = 10.0)
    groups, _ = _segment_loads_with_flag(metas; gap_k = gap_k)
    return groups
end

"""
    _segment_loads_with_flag(metas; gap_k = 10.0) -> (loads, flag)

`flag ∈ {:ok, :unimodal_fallback}`:
- `:ok`                — bimodal split found, or trivially one exposure
- `:unimodal_fallback` — no gap exceeds the threshold; entire directory is one load
"""
function _segment_loads_with_flag(metas::Vector{ExposureMeta}; gap_k::Float64 = 10.0)
    isempty(metas) && return (Vector{ExposureMeta}[], :ok)

    # Sort by timestamp; push exposures with missing timestamp to the end.
    sorted = sort(metas; by = m -> begin
        ts = m.prp !== nothing ? m.prp.timestamp : missing
        ismissing(ts) ? DateTime(9999) : ts
    end)

    length(sorted) == 1 && return ([[sorted[1]]], :ok)

    # Compute consecutive time gaps in seconds.
    timestamps = [begin
        ts = m.prp !== nothing ? m.prp.timestamp : missing
        ismissing(ts) ? missing : ts
    end for m in sorted]

    gaps = Float64[]
    for i in 2:length(timestamps)
        if !ismissing(timestamps[i]) && !ismissing(timestamps[i-1])
            push!(gaps, Float64(Dates.value(timestamps[i] - timestamps[i-1])) / 1000.0)  # ms → s
        end
    end

    if isempty(gaps)
        return ([sorted], :unimodal_fallback)
    end

    med = let v = sort(gaps)
        n = length(v)
        iseven(n) ? (v[n÷2] + v[n÷2+1]) / 2.0 : v[(n+1)÷2]
    end
    threshold = med * gap_k
    flag = :unimodal_fallback  # default; set to :ok when we actually split

    loads = Vector{ExposureMeta}[[sorted[1]]]
    for i in 2:length(sorted)
        ts_prev = sorted[i-1].prp !== nothing ? sorted[i-1].prp.timestamp : missing
        ts_curr = sorted[i].prp   !== nothing ? sorted[i].prp.timestamp   : missing
        gap_s = (!ismissing(ts_prev) && !ismissing(ts_curr)) ?
            Float64(Dates.value(ts_curr - ts_prev)) / 1000.0 : 0.0
        if gap_s > threshold
            push!(loads, ExposureMeta[])
            flag = :ok
        end
        push!(last(loads), sorted[i])
    end

    return (loads, flag)
end

# ---------------------------------------------------------------------------
# Stepping-axis detection + slot clustering (spec §5, step 1)
# ---------------------------------------------------------------------------

"""
    _median_inline(xs) -> Float64

Exact median of a non-empty `Vector{Float64}`. Inlined to avoid a `Statistics`
dependency (mirrors the inlined median in `_segment_loads_with_flag`).
"""
function _median_inline(xs::Vector{Float64})::Float64
    v = sort(xs)
    n = length(v)
    iseven(n) ? (v[n÷2] + v[n÷2+1]) / 2.0 : v[(n+1)÷2]
end

"""
    _cluster_slots(load_metas; slot_k = 5.0) -> Vector{Vector{ExposureMeta}}

Group the exposures in one load into slots (one slot = one sample position).

Algorithm (spec §5):
1. Extract horizontal positions; if all are missing → one slot (fallback).
2. Compute consecutive-pair position deltas; the within-burst jitter tolerance
   is `slot_k × median(|consecutive deltas|)`.
3. A position gap exceeding that tolerance starts a new slot.
4. Median-near-zero fallback (spec §5 "fall back to the absolute-position gap
   learned from the multi-frame bursts"): when the *median* consecutive delta is ~0
   — i.e. most consecutive frames sit at the *same* position (multi-frame bursts /
   within-slot revisits) — the plain `median × slot_k` tolerance degenerates to 0
   and would split on every step. Instead learn the tolerance from the non-zero
   deltas, which in this regime ARE the burst→burst (slot-to-slot) jumps. The
   tolerance must sit *between* the (≈0) within-burst jitter and the slot spacing,
   so it is `median(nonzero deltas) / slot_k` — a fraction of the learned slot
   spacing, below the jumps (so they split) yet above the jitter (so bursts don't).
   If there are no non-zero deltas at all (every frame at one position),
   tolerance = `Inf` → one slot.

`slot_k = 5.0` is a validated default. In the normal branch (non-zero median):
real data has 0.30 mm within-burst jitter vs 3.95 mm slot spacing → ratio 13, so a
5× multiplier on the jitter (1.5 mm) sits well inside the gap. In the fallback
branch the same `slot_k` divides the learned spacing (≈4 mm / 5 = 0.8 mm), which
likewise sits well below the jumps and above the (zero) jitter.

NOTE (deviation from the Phase-B plan / Task-6 brief draft): both drafts wrote the
fallback tolerance as `median(nonzero) × slot_k`. That is an authoring error — a
×5 tolerance (≈20 mm) exceeds the ~4 mm burst jumps it must split on, so it
collapses every burst-fixture load to a single slot and FAILS the brief's own
"3 bursts → 3 slots" assertion. Spec §5 requires the tolerance to fall *between*
jitter and slot spacing; `/slot_k` is the only reading consistent with both the
spec and the brief's behavioral test, so the operator is corrected to `/` here.

KNOWN LIMITATION (deferred — see spec §5 "single-frame acquisitions … else flag"):
a load of pure *single-frame* acquisitions visited round-robin (each slot exactly
one frame, no multi-frame burst anywhere) has every consecutive delta equal to a
slot-to-slot jump, so the jitter population and the spacing population are
indistinguishable — there is nothing to "learn the gap from". This case is NOT
handled here (it would under-split to one slot); the spec's prescribed behavior is
to flag for human review. Wiring that flag is deferred to Phase D's grouping-review
UI; for Phase B it falls through to a single slot.
"""
function _cluster_slots(
    load_metas::Vector{ExposureMeta};
    slot_k::Float64 = 5.0,
)::Vector{Vector{ExposureMeta}}
    isempty(load_metas) && return Vector{ExposureMeta}[]
    length(load_metas) == 1 && return [[load_metas[1]]]

    hpos = Union{Float64, Missing}[
        (m.prp !== nothing && !ismissing(m.prp.horizontal_position_mm)) ?
            m.prp.horizontal_position_mm : missing
        for m in load_metas
    ]

    if all(ismissing, hpos)
        # No position axis available: one slot, flag for review
        return [load_metas]
    end

    # Consecutive position deltas (absolute)
    deltas = Float64[]
    for i in 2:length(hpos)
        if !ismissing(hpos[i]) && !ismissing(hpos[i-1])
            push!(deltas, abs(hpos[i] - hpos[i-1]))
        end
    end

    if isempty(deltas)
        return [load_metas]
    end

    med_delta = _median_inline(deltas)
    # Median-near-zero fallback: when most consecutive frames sit at the same position
    # (multi-frame bursts), the median delta is ~0 and `median × slot_k` would be a
    # degenerate (zero) tolerance. Learn the slot-spacing tolerance from the non-zero
    # deltas (the burst→burst jumps) instead; Inf when there are no jumps at all.
    local threshold::Float64
    if med_delta < 1e-6
        nonzero = filter(d -> d > 1e-6, deltas)
        threshold = isempty(nonzero) ? Inf : _median_inline(nonzero) / slot_k
    else
        threshold = med_delta * slot_k
    end

    slots = Vector{ExposureMeta}[[load_metas[1]]]
    for i in 2:length(load_metas)
        h_prev = hpos[i-1]
        h_curr = hpos[i]
        gap = (!ismissing(h_prev) && !ismissing(h_curr)) ? abs(h_curr - h_prev) : 0.0
        if gap > threshold
            push!(slots, ExposureMeta[])
        end
        push!(last(slots), load_metas[i])
    end
    return slots
end

function _minimal_scan_config(
    tif_pattern::String,
    prp_pattern::String,
    dat_pattern::String,
    data_dir::String,
    analysis_dir::String,
)::ExperimentConfig
    # Positional constructor — ExperimentConfig has no @kwdef / keyword form.
    ExperimentConfig(
        # [experiment]
        "", "", "",
        # [beamline]
        nothing, nothing, "A-1", nothing, nothing, nothing,
        # [manifest]
        ",", 0, 0, 1, 1, 1, 1, 1, 1,
        # [layout]
        data_dir, analysis_dir, "simple",
        # [files]  integration_pattern=dat, image_pattern=tif
        dat_pattern, tif_pattern,
    )
end
