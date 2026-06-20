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

"""
    scan_directory(data_dir, analysis_dir;
                   tif_pattern = "{name}.tif",
                   prp_pattern = "{name}.prp",
                   dat_pattern = "{name}.dat") -> Vector{ExposureMeta}

Enumerate every TIFF found in `data_dir`, then pair each with its PRP and .dat
sidecar (the file the pattern names with `{name}` replaced by the stem). Returns
one `ExposureMeta` per TIFF stem, sorted by stem.

TIF enumeration scans the directory directly and strips the TIF suffix to derive
stems, rather than using `resolve_files` with an empty prefix (which
`_matches_prefix_with_boundary` rejects for alphanumeric-leading names — the
common case).
"""
function scan_directory(
    data_dir::AbstractString,
    analysis_dir::AbstractString;
    tif_pattern::String = "{name}.tif",
    prp_pattern::String = "{name}.prp",
    dat_pattern::String = "{name}.dat",
)::Vector{ExposureMeta}
    # Resolve a sidecar by substituting the known stem into the pattern: the file
    # is exactly `joinpath(base, pattern-with-{name}→stem)` if it exists.
    sidecar(base, stem, pattern) =
        (p = joinpath(base, replace(pattern, "{name}" => stem)); isfile(p) ? p : nothing)

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
        prp_path  = sidecar(data_dir, stem, prp_pattern)
        dat_path  = sidecar(analysis_dir, stem, dat_pattern)
        prp_parsed = prp_path !== nothing ? parse_prp(prp_path) : nothing
        push!(metas, ExposureMeta(stem, tif_path, prp_path, dat_path, prp_parsed))
    end
    return metas
end

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

# ---------------------------------------------------------------------------
# Grouped structs returned by group_into_samples (spec §5)
# ---------------------------------------------------------------------------

struct GroupedExposure
    stem       ::String
    tif_path   ::Union{String, Nothing}
    prp_path   ::Union{String, Nothing}
    dat_path   ::Union{String, Nothing}
    timestamp  ::Union{DateTime, Missing}
    exposure_time_s        ::Union{Float64, Missing}
    horizontal_position_mm ::Union{Float64, Missing}
    scan_id    ::Union{Int, Missing}   # parsed from the `_S<digits>_` filename token
    frame_no   ::Union{Int, Missing}   # parsed from the trailing `_<digits>` frame index
end

struct GroupedSample
    name           ::String
    name_source    ::String   # always "auto" from this function
    slot_index     ::Int
    grouping_source::String   # "auto_position" or "auto_time"
    exposures      ::Vector{GroupedExposure}
end

struct GroupedLoad
    load_index  ::Int
    frame_count ::Int
    start_time  ::Union{DateTime, Missing}
    end_time    ::Union{DateTime, Missing}
    samples     ::Vector{GroupedSample}
    flag        ::Symbol   # :ok | :unimodal_fallback | :single_exposure
end

struct GroupingResult
    loads        ::Vector{GroupedLoad}
    discrepancies::Vector{String}  # human-readable grouping anomalies
end

# ---------------------------------------------------------------------------
# Auto-naming helper (spec §5 naming rule)
# ---------------------------------------------------------------------------

"""
    _auto_name(label_hint, load_index, slot_index) -> String

Produce a sample name per spec §5: `<label> (SNNPMM)` where S=load index,
P=slot index (both 1-based, zero-padded to two digits). `label_hint` is derived
from the filename token when parseable; an empty hint yields just the coordinate
anchor `(SNNPMM)`.
"""
function _auto_name(label_hint::AbstractString, load_index::Int, slot_index::Int)
    coord = "S$(lpad(load_index, 2, '0'))P$(lpad(slot_index, 2, '0'))"
    isempty(strip(label_hint)) && return "($(coord))"
    return "$(strip(label_hint)) ($(coord))"
end

"""
    _label_from_stem(stem) -> String

Extract the human label token from a filename stem like `HA_85_422_S2404_0_001`.
Heuristic: the portion BEFORE the scan-id token `_S\\d+_`. Returns the stem
unchanged if no such token is found.
"""
function _label_from_stem(stem::AbstractString)
    m = match(r"^(.+?)_S\d+_", stem)
    m === nothing && return String(stem)
    return String(m.captures[1])
end

"""
    _parse_scan_frame(stem) -> (scan_id::Union{Int,Missing}, frame_no::Union{Int,Missing})

Extract the SSRL scan id and frame index from a filename stem such as
`HA_85_422_S2404_0_001` (scan_id=2404, frame_no=1). The scan id is the integer
in the `_S<digits>_` token; the frame number is the trailing `_<digits>` group
at the very end of the stem. Either is `missing` when its token is absent
(confirmed against the real SSRL 2026-04 naming convention).
"""
function _parse_scan_frame(stem::AbstractString)
    m_scan  = match(r"_S(\d+)_", stem)
    m_frame = match(r"_(\d+)$", stem)
    scan_id  = m_scan  === nothing ? missing : parse(Int, m_scan.captures[1])
    frame_no = m_frame === nothing ? missing : parse(Int, m_frame.captures[1])
    return (scan_id, frame_no)
end

# ---------------------------------------------------------------------------
# group_into_samples (spec §5)
# ---------------------------------------------------------------------------

"""
    group_into_samples(metas::Vector{ExposureMeta}) -> GroupingResult

Run the full §5 backbone on the flat list of per-exposure metadata:
  1. Time-gap segment into loads (`_segment_loads_with_flag`).
  2. Within each load, cluster by position (`_cluster_slots`).
  3. Build the `GroupedLoad`/`GroupedSample`/`GroupedExposure` tree with
     auto-names (`<label> (SNNPMM)`).

The segmentation fallback (`:unimodal_fallback`) is surfaced both as a
human-readable discrepancy string and onto each affected load's `flag`.
"""
# ---------------------------------------------------------------------------
# Read-time sample-flag derivation (spec §8.8 / §9.1)
#
# PURE: takes the get_loads_rollup rows (Phase D), returns a Dict keyed by
# sample_id. NO DB, never a stored column. Phase D's get_loads_rollup calls this
# over the rows it reads and suppresses any flag with a non-undone
# grouping_flag_dismissed event before serializing (spec §9.2/§9.3).
# ---------------------------------------------------------------------------

"A cross-load merge suggestion: this sample's filename label recurs in another load."
struct MergeFlag
    merge_with_sample_id ::Int
    merge_with_label     ::String
end

"An intra-sample split suggestion: the sample spans a position jump beyond local jitter."
struct SplitFlag
    split_at_index ::Int        # 1-based exposure index where the jump occurs
    jump_from      ::Float64    # position just before the jump
    jump_to        ::Float64    # position at the jump
end

"""One flag per sample, serialized by Phase D to the §8.8 `GroupingFlag` JSON union."""
const GroupingFlag = Union{MergeFlag, SplitFlag}

"""
    derive_sample_flags(load_rows; min_slot_separation_mm = 0.5) -> Dict{Int, GroupingFlag}

PURE read-time derivation of per-sample merge/split suggestions over the
`get_loads_rollup` rows (spec §8.8). No DB access; returns a Dict keyed by
`sample_id`, with a sample **absent** from the Dict meaning "no flag" (the
contract's JSON `null`).

`load_rows` is the `get_loads_rollup` shape (verified 2026-06-18):
`Vector` of loads `(load_id, load_index, …, samples)`; each sample
`(sample_id, name, slot_index, …, exposures)`; each exposure
`(id, filename, horizontal_position, timestamp)` (the §8.8 leaf).

Two suggestion kinds (a sample gets at most one; **split wins** if both apply):

1. **Split** — within one sample, a consecutive `horizontal_position` gap
   exceeding `min_slot_separation_mm` marks a genuine slot boundary (split
   candidate). Gaps at or below it are stage jitter (no split).

   Physical rationale: SAXS samples sit in quartz capillaries (0.5–2 mm
   diameter) and the X-ray beam passes through exactly one at a time, so
   adjacent slots are physically ≥ ~1 mm apart. Stage repeatability when
   returning to the *same* slot is ~0.1–0.3 mm in real SSRL data. A gap of
   0.5 mm therefore sits ~5× above stage jitter and ~2× below the minimum
   slot pitch — a physical floor, not a fitted constant. It errs toward
   flagging (false flags are cheap; these are human-review suggestions,
   nothing auto-changes).

2. **Merge** — a sample's filename label (via `_label_from_stem` over its
   exposures' `filename`s) recurs as the label of a sample in *another* load. The
   grouper never auto-merges cross-load (spec §5), so this is surfaced as a
   suggestion pointing at the other sample (`merge_with_sample_id` +
   `merge_with_label`). When a label recurs in more than two loads, each flagged
   sample points at the *first other* sample sharing that label (lowest
   `sample_id`); the UI walks the chain.
"""
function derive_sample_flags(load_rows;
        min_slot_separation_mm::Float64 = 0.5)::Dict{Int, GroupingFlag}
    flags = Dict{Int, GroupingFlag}()

    # =====================================================================
    # 1. SPLIT suggestions — per sample: any consecutive position gap above
    #    min_slot_separation_mm is a genuine slot boundary; anything at or
    #    below is stage return jitter. Single absolute physical threshold —
    #    no median computation, no regime switching.
    #
    #    Physical basis: SAXS capillary slots are ≥ ~1 mm apart; stage
    #    repeatability (returning to the same slot) is ~0.1–0.3 mm. The 0.5 mm
    #    default sits ~5× above jitter and ~2× below the minimum slot pitch.
    # =====================================================================
    for ld in load_rows
        for sm in ld.samples
            xs = sm.exposures
            length(xs) < 2 && continue
            for i in 2:length(xs)
                p_prev = xs[i-1].horizontal_position
                p_curr = xs[i].horizontal_position
                # Skip unless BOTH positions are real numbers. `isa Real` catches
                # `nothing` AND SQLite NULL (`missing`) in one guard — a bare
                # `=== nothing` would let `missing` through and Float64(missing) throws.
                (p_prev isa Real && p_curr isa Real) || continue
                if abs(Float64(p_curr) - Float64(p_prev)) > min_slot_separation_mm
                    flags[Int(sm.sample_id)] = SplitFlag(
                        i, Float64(p_prev), Float64(p_curr))
                    break  # first jump only; the UI resolves one split at a time
                end
            end
        end
    end

    # =====================================================================
    # 2. MERGE suggestions — a sample's filename label recurs in another load.
    #    Compute each sample's label once, group sample_ids by (label), and flag
    #    every member of a multi-load group. Split takes precedence (skip if the
    #    sample already has a split flag).
    # =====================================================================
    # label -> Vector of (sample_id, load_id), in stable iteration order.
    by_label = Dict{String, Vector{Tuple{Int, Int}}}()
    for ld in load_rows
        for sm in ld.samples
            isempty(sm.exposures) && continue
            label = _label_from_stem(first(sm.exposures).filename)
            push!(get!(by_label, label, Tuple{Int, Int}[]),
                  (Int(sm.sample_id), Int(ld.load_id)))
        end
    end

    for (label, members) in by_label
        # Recurrence requires the label to appear in >1 distinct load.
        distinct_loads = unique(last.(members))
        length(distinct_loads) < 2 && continue
        for (sid, lid) in members
            haskey(flags, sid) && continue  # split wins
            # Point at the first OTHER sample with this label in a different load
            # (lowest sample_id for determinism).
            others = sort([s for (s, l) in members if s != sid])
            isempty(others) && continue
            flags[sid] = MergeFlag(first(others), label)
        end
    end

    return flags
end

function group_into_samples(metas::Vector{ExposureMeta})::GroupingResult
    isempty(metas) && return GroupingResult(GroupedLoad[], String[])

    load_groups, seg_flag = _segment_loads_with_flag(metas)
    discrepancies = String[]

    # Surface the segmentation fallback durably (spec §5: "unimodal fallback → flag
    # for review"). The flag is a property of the whole segmentation pass, so it is
    # both (a) recorded as a human-readable discrepancy and (b) carried onto the load
    # row(s) below via `seg_flag`.
    if seg_flag == :unimodal_fallback
        push!(discrepancies,
            "load segmentation: no clear bimodal time-gap separation; treated the " *
            "entire directory as a single load — review load boundaries")
    end

    grouped_loads = GroupedLoad[]
    for (li, load_metas) in enumerate(load_groups)
        slot_groups = _cluster_slots(load_metas)
        grouped_samples = GroupedSample[]
        for (si, slot_metas) in enumerate(slot_groups)
            # Label hint drawn from the slot's first stem.
            label_hint = _label_from_stem(first(slot_metas).stem)

            exps = [begin
                sid, fno = _parse_scan_frame(m.stem)
                ts  = (m.prp !== nothing && !ismissing(m.prp.timestamp)) ? m.prp.timestamp : missing
                ets = m.prp !== nothing ? m.prp.exposure_time_s : missing
                hpm = m.prp !== nothing ? m.prp.horizontal_position_mm : missing
                GroupedExposure(
                    m.stem,
                    m.tif_path, m.prp_path, m.dat_path,
                    ts, ets, hpm,
                    sid, fno,
                )
            end for m in slot_metas]

            push!(grouped_samples, GroupedSample(
                _auto_name(label_hint, li, si),
                "auto",
                si,
                "auto_position",
                exps,
            ))
        end

        timestamps = DateTime[]
        for m in load_metas
            ts = m.prp !== nothing ? m.prp.timestamp : missing
            ismissing(ts) || push!(timestamps, ts)
        end
        start_time = isempty(timestamps) ? missing : minimum(timestamps)
        end_time   = isempty(timestamps) ? missing : maximum(timestamps)

        push!(grouped_loads, GroupedLoad(
            li,
            length(load_metas),
            start_time,
            end_time,
            grouped_samples,
            seg_flag,   # :ok | :unimodal_fallback | :single_exposure — from segmentation (spec §5)
        ))
    end

    return GroupingResult(grouped_loads, discrepancies)
end
