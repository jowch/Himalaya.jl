# Ingestion Redesign — Phase B: Ingest Core Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the four new backend modules — `prp.jl`, `geometry.jl`, `grouping.jl`, and the `scan_and_group!` orchestrator in a new `ingest.jl` — that together replace the manifest-driven `_reingest_inner!` path with a fully automatic "point at a directory → scan, read PRP/setup files, group into loads and samples, analyze" pipeline. After this phase the CLI `init`/`reingest` commands can call the same core the HTTP API will call in Phase C.

**Architecture:** Four layered modules, each independently testable:
1. `prp.jl` — pure parser; no DB, no I/O beyond the single file.
2. `geometry.jl` — pure aggregator over parsed PRPs + setup files; returns derived geometry + discrepancy list; no DB.
3. `grouping.jl` — pure grouper over per-exposure metadata; returns load/sample/exposure structs; no DB. Reuses `resolve_files` from `config.jl` for filesystem discovery.
4. `ingest.jl` — transactional orchestrator; calls `prp.jl` → `geometry.jl` → `grouping.jl` → DB writers from Phase A (`create_experiment!`, `create_sample!`, `create_exposure!`, `analyze_exposure!`); clones `_reingest_inner!`'s insert-only discipline + `_DB_WRITE_LOCK` from `server.jl:29`.

**Tech Stack:** Julia, stdlib `Test`, `Dates`, `Statistics`. No new package dependencies. Backend package at `packages/HimalayaUI/`.

**Spec:** `docs/superpowers/specs/2026-06-15-ingestion-redesign-design.md` §5 (grouping algorithm), §6 (geometry derivation), §9.1 (shared grouping core). Read those sections before starting; §5 in particular is the source of truth for all threshold decisions.

**Source of truth for current code:** the build-kit anchors in this plan were line-verified 2026-06-18, but line numbers drift — confirm each anchor with a quick `grep`/read before editing.

---

## File Structure

| File | Responsibility | This plan |
|---|---|---|
| `packages/HimalayaUI/src/prp.jl` | `parse_prp(path) → NamedTuple`: key=value PRP parser | CREATE |
| `packages/HimalayaUI/src/geometry.jl` | `derive_geometry(prp_paths, setup_files) → (geometry, discrepancies)` | CREATE |
| `packages/HimalayaUI/src/grouping.jl` | `scan_directory`, `group_into_samples` | CREATE |
| `packages/HimalayaUI/src/ingest.jl` | `scan_and_group!(db, exp_id, dir; additive)` orchestrator + `cheap_change_check(db, exp_id, dir)` (Phase-C contract) | CREATE |
| `packages/HimalayaUI/src/HimalayaUI.jl` | module includes | MODIFY (add four new includes) |
| `packages/HimalayaUI/test/test_ingestion_core.jl` | new standalone test file for this phase | CREATE |
| `packages/HimalayaUI/test/runtests.jl` | test registry | MODIFY (add include) |

**Out of scope (later plans):** the Phase-A schema migrations (Phase A); the HTTP API routes `POST /api/experiments`, `POST /api/experiments/{id}/scan`, `GET /api/experiments/{id}/loads`, and the `broadcast_progress!` SSE helper (Phase C); the structural-edit event kinds `exposure_moved`/`sample_renamed`/`sample_created`/`sample_split`/`grouping_flag_dismissed` (Phase D; no `sample_merged`); the frontend `/experiments/:id` shell and all new React components (Phase E). The CLI `init`/`reingest` rewrite to call `scan_and_group!` is the last step of this phase (Task 10), not Phase C, because it reuses the same core and can be done without the HTTP layer.

---

## Conventions for every task

- **Run a single test file** during TDD from the repo root:
  `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
  (Standalone: the file ends with its own `@testset`s and runs directly without the full suite.)
- Each test builds **in-memory fixtures** (synthetic PRP/setup text written to `mktempdir()`) so tests run offline, never touching the real SSRL volume at `/Volumes/data`.
- Tests that need a DB use `HimalayaUI.open_db(joinpath(mktempdir(), "h.db"))` — the standard pattern confirmed across the existing test suite.
- **Commit after each task** once its test passes.
- New modules are added to `HimalayaUI.jl` immediately when created (the module won't compile without the include).

---

## Task 0: Test harness + helpers

**Files:**
- Create: `packages/HimalayaUI/test/test_ingestion_core.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl`

- [ ] **Step 1: Create the standalone test file with fixture helpers**

```julia
# packages/HimalayaUI/test/test_ingestion_core.jl
using Test
using SQLite, DBInterface, Tables
using Dates
using HimalayaUI

# ---------------------------------------------------------------------------
# Shared fixture helpers
# ---------------------------------------------------------------------------

"""
Write a synthetic PRP file to `path` with the given field overrides.
Defaults match the real HA_85_422_S2404_0_001.prp from the SSRL 2026-04 run.

Note: scan id and frame number are NOT PRP fields — they live in the *filename
stem* (`_S<digits>_…_<frame>`) and are recovered by `_parse_scan_frame` (Task 7),
so `write_prp` takes no `scan_id`/`frame_no` kwargs. Tests exercise scan/frame
parsing by naming the fixture file accordingly (e.g. `HA_85_422_S2404_0_001.prp`).
"""
function write_prp(path::AbstractString;
        timestamp = "26 Apr 2026 18:20:27",
        beam_energy_ev = 9000.027604502573,
        pipe_length_mm = 1700,
        detector = "Pilatus 1M",
        exposure_time = 15,
        horizontal_position_mm = 58.9)
    open(path, "w") do io
        println(io, "Image file name: $(basename(path)[1:end-4]).tif")
        println(io, "Time this file was written: $timestamp")
        println(io, "Exposure time=$exposure_time")
        println(io, "Counting time=$(exposure_time + 0.8)")
        println(io, "Beam energy=$beam_energy_ev eV")
        println(io, "Pipe length=$pipe_length_mm mm")
        println(io, "Scan motor=motor")
        println(io, "Scan range=delta mm")
        println(io, "Horizontal position=$horizontal_position_mm mm")
        println(io, "Vertical position=26.55 mm")
        println(io, "dispx position=9.8 mm")
        println(io, "dispy position=35.9 mm")
        println(io, "Detector=$detector")
        println(io, "Detector mode=Pilatus; 0 for normal and 1 for dezingered")
        println(io, "Phi position=157")
        println(io, "Attenuator=0 um Al")
    end
end

"""
Write a synthetic setup_info file to `path`.
Defaults match the real setup_info_20260425_181705.txt from SSRL 2026-04.
"""
function write_setup_info(path::AbstractString;
        beam_center_x = 421.409,
        beam_center_y = 836.946,
        mean_distance_mm = 1809.5,
        ring_distances = [1806.0, 1811.1, 1811.5])
    open(path, "w") do io
        println(io, "File used for calibration:")
        println(io, "../data/agbe_1p7m_S1958_0_001.tif")
        println(io, "")
        println(io, "Beam center is at ($beam_center_x, $beam_center_y)")
        println(io, "")
        println(io, "Calculated detector distances are:")
        for (i, d) in enumerate(ring_distances)
            println(io, "Ring $i:\t pixel $(247.82 * i):\t $d mm")
        end
        println(io, "----------------------------------------------------------")
        println(io, "Mean distance:\t\t $mean_distance_mm mm")
    end
end

"Open a fresh DB with the full migrated schema."
function fresh_db()
    path = joinpath(mktempdir(), "h.db")
    HimalayaUI.open_db(path)
end

@testset "ingestion core (Phase B)" begin
    # task @testsets are appended below
end
```

- [ ] **Step 2: Register in runtests.jl**

Inside the `@testset "HimalayaUI"` block, add the include after the actual last entry (verified 2026-06-18: the final two includes are `test_migrate_toml.jl` then `test_spa_fallback.jl`) and before the closing `end`:

```julia
    include("test_ingestion_core.jl")
```

(Placement within the block is not load-bearing — anywhere inside the `@testset` works; appending after `test_spa_fallback.jl` keeps it last.)

- [ ] **Step 3: Run to confirm harness loads**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: PASS (empty `@testset` passes; proves imports resolve).

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/test/test_ingestion_core.jl packages/HimalayaUI/test/runtests.jl
git commit -m "test: scaffold ingestion core test harness"
```

---

## Task 1: `prp.jl` — `parse_prp`

Parse a single PRP file into a `NamedTuple`. The format is `Key=Value` lines (confirmed against the real `HA_85_422_S2404_0_001.prp`). Per spec §9.1: defensive — unparseable fields become `missing`. Return only the fields the system uses; drop motor slits, attenuator, and per-counter values (spec §4: "most PRP fields are dropped").

**Files:**
- Create: `packages/HimalayaUI/src/prp.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl` (add `include("prp.jl")`)
- Test: `packages/HimalayaUI/test/test_ingestion_core.jl`

- [ ] **Step 1: Write the failing test**

Append to the `@testset "ingestion core (Phase B)"` block:

```julia
@testset "parse_prp" begin
    dir = mktempdir()
    prp_path = joinpath(dir, "HA_85_422_S2404_0_001.prp")
    write_prp(prp_path;
        timestamp = "26 Apr 2026 18:20:27",
        beam_energy_ev = 9000.027604502573,
        pipe_length_mm = 1700,
        detector = "Pilatus 1M",
        exposure_time = 15.0,
        horizontal_position_mm = 58.9)
    prp = HimalayaUI.parse_prp(prp_path)

    # Energy: strip "eV", convert to keV
    @test prp.beam_energy_ev ≈ 9000.027604502573
    @test prp.energy_kev ≈ 9.000027604502573

    # Pipe length: strip "mm", store as metres
    @test prp.pipe_length_m ≈ 1.700

    # Detector model string (drives pitch lookup in geometry.jl)
    @test prp.detector == "Pilatus 1M"

    # Exposure time (seconds)
    @test prp.exposure_time_s ≈ 15.0

    # Horizontal position (mm)
    @test prp.horizontal_position_mm ≈ 58.9

    # Timestamp parses to a DateTime
    @test prp.timestamp isa DateTime
    @test year(prp.timestamp) == 2026
    @test month(prp.timestamp) == 4

    # Graceful missing: a truncated PRP returns missing for absent fields
    trunc_path = joinpath(dir, "trunc.prp")
    write(trunc_path, "Image file name: trunc.tif\nBeam energy=9000 eV\n")
    trunc_prp = HimalayaUI.parse_prp(trunc_path)
    @test trunc_prp.horizontal_position_mm === missing
    @test trunc_prp.detector === missing
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: FAIL (`UndefVarError: parse_prp not defined`).

- [ ] **Step 3: Create `prp.jl`**

```julia
# packages/HimalayaUI/src/prp.jl
#
# PRP file parser (§9.1 shared ingest core).
# Format: "Key=Value" lines with optional trailing unit tokens (mm, eV, etc.).
# Confirmed against the real SSRL 2026-04 PRP: HA_85_422_S2404_0_001.prp.
#
# Returned NamedTuple fields:
#   timestamp              :: Union{DateTime, Missing}   — "Time this file was written:"
#   beam_energy_ev         :: Union{Float64, Missing}    — raw eV value
#   energy_kev             :: Union{Float64, Missing}    — beam_energy_ev / 1000
#   pipe_length_m          :: Union{Float64, Missing}    — "Pipe length" converted m
#   detector               :: Union{String, Missing}     — e.g. "Pilatus 1M"
#   exposure_time_s        :: Union{Float64, Missing}    — "Exposure time"
#   horizontal_position_mm :: Union{Float64, Missing}    — "Horizontal position"
#
# Dropped fields (constant/noise per spec §4): Vertical position, dispx, dispy,
# detector_vert/horz, Phi, slit widths, attenuator, per-counter voltages, scan motor.

using Dates

"""
    parse_prp(path) -> NamedTuple

Parse one SSRL PRP file at `path` and return the fields the ingest system uses.
Unparseable or absent fields are `missing`. Never throws on malformed lines;
a `@warn` is emitted instead so batch ingestion can proceed.
"""
function parse_prp(path::AbstractString)
    isfile(path) || error("PRP file not found: $path")

    timestamp              = missing
    beam_energy_ev         = missing
    energy_kev             = missing
    pipe_length_m          = missing
    detector               = missing
    exposure_time_s        = missing
    horizontal_position_mm = missing

    for line in eachline(path)
        # Timestamp line has a different format: "Time this file was written: DD Mon YYYY HH:MM:SS"
        if startswith(line, "Time this file was written:")
            raw = strip(line[length("Time this file was written:") + 1:end])
            try
                timestamp = DateTime(raw, dateformat"dd u yyyy HH:MM:SS")
            catch
                @warn "parse_prp: could not parse timestamp" path line
            end
            continue
        end

        idx = findfirst('=', line)
        idx === nothing && continue
        key = strip(line[1:idx-1])
        val_raw = strip(line[idx+1:end])

        # Strip trailing unit token (first non-numeric, non-dot, non-sign char onward).
        # E.g. "9000.03 eV" → 9000.03; "1700 mm" → 1700; "15" → 15.
        val_str = let s = val_raw
            m = match(r"^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?", s)
            m === nothing ? nothing : m.match
        end

        if key == "Beam energy" && val_str !== nothing
            try
                beam_energy_ev = parse(Float64, val_str)
                energy_kev     = beam_energy_ev / 1000.0
            catch
                @warn "parse_prp: could not parse Beam energy" path val_raw
            end

        elseif key == "Pipe length" && val_str !== nothing
            try
                pipe_length_m = parse(Float64, val_str) / 1000.0
            catch
                @warn "parse_prp: could not parse Pipe length" path val_raw
            end

        elseif key == "Detector"
            # e.g. "Pilatus 1M" — string, no numeric conversion needed
            detector = val_raw

        elseif key == "Exposure time" && val_str !== nothing
            try
                exposure_time_s = parse(Float64, val_str)
            catch
                @warn "parse_prp: could not parse Exposure time" path val_raw
            end

        elseif key == "Horizontal position" && val_str !== nothing
            try
                horizontal_position_mm = parse(Float64, val_str)
            catch
                @warn "parse_prp: could not parse Horizontal position" path val_raw
            end
        end
    end

    return (
        timestamp              = timestamp,
        beam_energy_ev         = beam_energy_ev,
        energy_kev             = energy_kev,
        pipe_length_m          = pipe_length_m,
        detector               = detector,
        exposure_time_s        = exposure_time_s,
        horizontal_position_mm = horizontal_position_mm,
    )
end
```

- [ ] **Step 4: Add to HimalayaUI.jl**

Read `packages/HimalayaUI/src/HimalayaUI.jl`. Insert `include("prp.jl")` immediately before `include("config.jl")` (line 5) so prp.jl's Dates import is visible before config:

```julia
include("db.jl")
include("image.jl")
include("datfile.jl")
include("prp.jl")        # ← new
include("config.jl")
```

- [ ] **Step 5: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/prp.jl packages/HimalayaUI/src/HimalayaUI.jl packages/HimalayaUI/test/test_ingestion_core.jl
git commit -m "feat(prp): add parse_prp key=value PRP parser"
```

---

## Task 2: `geometry.jl` — detector pitch lookup + `parse_setup_info`

Two sub-components assembled here: the static detector→pitch table (spec §6: "a static table in `geometry.jl`") and the setup-file parser. The `derive_geometry` orchestrator that calls both comes in Task 3.

**Files:**
- Create: `packages/HimalayaUI/src/geometry.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_core.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "detector pitch lookup" begin
    # Known detector
    @test HimalayaUI.detector_pixel_size_um("Pilatus 1M") ≈ 172.0
    @test HimalayaUI.detector_pixel_size_um("Pilatus 2M") ≈ 172.0
    @test HimalayaUI.detector_pixel_size_um("Pilatus 300K") ≈ 172.0
    @test HimalayaUI.detector_pixel_size_um("Eiger 1M")  ≈ 75.0
    @test HimalayaUI.detector_pixel_size_um("Eiger 4M")  ≈ 75.0
    # Unknown detector → missing
    @test HimalayaUI.detector_pixel_size_um("Unknown XRD") === missing
end

@testset "parse_setup_info" begin
    dir = mktempdir()
    setup_path = joinpath(dir, "setup_info_20260425_181705.txt")
    write_setup_info(setup_path;
        beam_center_x = 421.409,
        beam_center_y = 836.946,
        mean_distance_mm = 1809.5)
    info = HimalayaUI.parse_setup_info(setup_path)
    @test info.beam_center_x ≈ 421.409
    @test info.beam_center_y ≈ 836.946
    @test info.mean_distance_m ≈ 1.8095

    # Missing Mean distance line → missing
    bad_path = joinpath(dir, "setup_bad.txt")
    write(bad_path, "Beam center is at (100.0, 200.0)\n")
    bad_info = HimalayaUI.parse_setup_info(bad_path)
    @test bad_info.beam_center_x ≈ 100.0
    @test bad_info.mean_distance_m === missing

    # Completely empty → all missing
    empty_path = joinpath(dir, "setup_empty.txt")
    write(empty_path, "")
    empty_info = HimalayaUI.parse_setup_info(empty_path)
    @test empty_info.beam_center_x === missing
    @test empty_info.mean_distance_m === missing
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: FAIL.

- [ ] **Step 3: Create `geometry.jl` with pitch lookup and setup parser**

```julia
# packages/HimalayaUI/src/geometry.jl
#
# Geometry derivation (spec §6):
#   1. detector_pixel_size_um(model) — static lookup table.
#   2. parse_setup_info(path) — parse one analysis/setup_info_*.txt file.
#   3. derive_geometry(prp_paths, setup_files) — orchestrate (Task 3).
#
# Setup-file format confirmed against real SSRL data (2026-06-18):
#   analysis/setup_info_20260425_181705.txt (11 lines):
#     "Beam center is at (421.409, 836.946)"
#     "Mean distance:         1809.5 mm"
#
# Geometry authority (spec §6, confirmed Jonathan 2026-06-18):
#   flight_path_m ← setup file Mean distance (AgBe calibration) when present,
#                   else PRP Pipe length.  Source tagged accordingly.

# ---------------------------------------------------------------------------
# Detector → pixel pitch lookup
# ---------------------------------------------------------------------------

"""Detector model string → pixel pitch in µm, or `missing` for unknown models."""
const _DETECTOR_PITCH_UM = Dict{String, Float64}(
    # Pilatus family (Dectris): all 172 µm pixels
    "Pilatus 1M"   => 172.0,
    "Pilatus 2M"   => 172.0,
    "Pilatus 6M"   => 172.0,
    "Pilatus 300K" => 172.0,
    "Pilatus 200K" => 172.0,
    # Eiger family (Dectris): all 75 µm pixels
    "Eiger 1M"     => 75.0,
    "Eiger 4M"     => 75.0,
    "Eiger 9M"     => 75.0,
    "Eiger 16M"    => 75.0,
)

"""
    detector_pixel_size_um(model) -> Union{Float64, Missing}

Look up the pixel pitch (µm) for a detector model string.
Returns `missing` for unknown models; the caller should flag for manual entry
rather than guessing (spec §6).
"""
function detector_pixel_size_um(model::AbstractString)::Union{Float64, Missing}
    # Exact match first
    haskey(_DETECTOR_PITCH_UM, model) && return _DETECTOR_PITCH_UM[model]
    # Prefix-based match for strings like "Pilatus 1M (dezingered)"
    for (k, v) in _DETECTOR_PITCH_UM
        startswith(model, k) && return v
    end
    return missing
end

# ---------------------------------------------------------------------------
# setup_info_*.txt parser
# ---------------------------------------------------------------------------

"""
    parse_setup_info(path) -> NamedTuple

Parse one `analysis/setup_info_<YYYYMMDD_HHMMSS>.txt` file and extract the
beam center and calibrated detector distance. Returns `missing` for absent fields.

Confirmed format (SSRL 2026-04, setup_info_20260425_181705.txt):
  "Beam center is at (421.409, 836.946)"
  "Mean distance:         1809.5 mm"
"""
function parse_setup_info(path::AbstractString)
    isfile(path) || error("setup_info file not found: $path")

    beam_center_x    = missing
    beam_center_y    = missing
    mean_distance_m  = missing

    for line in eachline(path)
        # "Beam center is at (X, Y)"
        m_bc = match(r"Beam center is at \(\s*([0-9.]+)\s*,\s*([0-9.]+)\s*\)", line)
        if m_bc !== nothing
            try
                beam_center_x = parse(Float64, m_bc.captures[1])
                beam_center_y = parse(Float64, m_bc.captures[2])
            catch
                @warn "parse_setup_info: could not parse beam center" path line
            end
            continue
        end

        # "Mean distance:         1809.5 mm"
        m_md = match(r"Mean distance:\s*([0-9.]+)\s*mm", line)
        if m_md !== nothing
            try
                mean_distance_m = parse(Float64, m_md.captures[1]) / 1000.0
            catch
                @warn "parse_setup_info: could not parse Mean distance" path line
            end
            continue
        end
    end

    return (
        beam_center_x   = beam_center_x,
        beam_center_y   = beam_center_y,
        mean_distance_m = mean_distance_m,
    )
end
```

- [ ] **Step 4: Add to HimalayaUI.jl**

Insert `include("geometry.jl")` immediately after `include("prp.jl")`:

```julia
include("prp.jl")
include("geometry.jl")   # ← new
include("config.jl")
```

- [ ] **Step 5: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/geometry.jl packages/HimalayaUI/src/HimalayaUI.jl packages/HimalayaUI/test/test_ingestion_core.jl
git commit -m "feat(geometry): detector pitch lookup + setup_info parser"
```

---

## Task 3: `geometry.jl` — `derive_geometry` orchestrator

Combines parsed PRPs + setup files into a single geometry result with per-field source tags and a discrepancy list (spec §6). Multi-setup variance detection is included: if PRP beam-energy or pipe-length varies across the directory's PRPs, that is surfaced as a discrepancy.

**Files:**
- Modify: `packages/HimalayaUI/src/geometry.jl` (append `derive_geometry`)
- Test: `packages/HimalayaUI/test/test_ingestion_core.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "derive_geometry" begin
    dir = mktempdir()
    data_dir = joinpath(dir, "data")
    mkpath(data_dir)
    analysis_dir = joinpath(dir, "analysis")
    mkpath(analysis_dir)

    # Two PRP files with consistent geometry
    for (name, hpos) in [("HA_001", 58.9), ("HA_002", 63.1)]
        write_prp(joinpath(data_dir, "$name.prp");
            beam_energy_ev = 9000.027604502573,
            pipe_length_mm = 1700,
            detector = "Pilatus 1M",
            exposure_time = 15.0,
            horizontal_position_mm = hpos)
    end

    # One setup file with calibrated distance
    write_setup_info(joinpath(analysis_dir, "setup_info_20260425_181705.txt");
        beam_center_x = 421.409, beam_center_y = 836.946, mean_distance_mm = 1809.5)

    prp_paths   = [joinpath(data_dir, "HA_001.prp"), joinpath(data_dir, "HA_002.prp")]
    setup_files = [joinpath(analysis_dir, "setup_info_20260425_181705.txt")]

    geo, discrepancies = HimalayaUI.derive_geometry(prp_paths, setup_files)

    # Energy from PRP
    @test geo.energy_kev ≈ 9.000027604502573
    @test geo.energy_kev_source == "prp"

    # flight_path_m from calibrated setup file (authority per spec §6)
    @test geo.flight_path_m ≈ 1.8095
    @test geo.flight_path_m_source == "setup"

    # Beam center from setup file
    @test geo.beam_center_x ≈ 421.409
    @test geo.beam_center_y ≈ 836.946
    @test geo.beam_center_x_source == "setup"

    # Pixel size from detector→pitch lookup
    @test geo.pixel_size_um ≈ 172.0
    @test geo.pixel_size_um_source == "prp"

    # No discrepancies when fields are consistent
    @test isempty(discrepancies)

    # Discrepancy when PRP beam energies disagree
    write_prp(joinpath(data_dir, "HA_003.prp");
        beam_energy_ev = 8900.0,    # different!
        pipe_length_mm = 1700, detector = "Pilatus 1M",
        horizontal_position_mm = 68.0)
    geo2, disc2 = HimalayaUI.derive_geometry(
        vcat(prp_paths, [joinpath(data_dir, "HA_003.prp")]), setup_files)
    @test any(d -> d.field == "beam_energy_ev", disc2)

    # No setup file → fall back to PRP pipe length
    geo3, _ = HimalayaUI.derive_geometry(prp_paths, String[])
    @test geo3.flight_path_m ≈ 1.700
    @test geo3.flight_path_m_source == "prp"

    # Unknown detector → pixel_size_um missing, discrepancy flagged
    write_prp(joinpath(data_dir, "HA_004.prp");
        beam_energy_ev = 9000.0, pipe_length_mm = 1700,
        detector = "SuperDetector X",
        horizontal_position_mm = 72.0)
    geo4, disc4 = HimalayaUI.derive_geometry([joinpath(data_dir, "HA_004.prp")], String[])
    @test geo4.pixel_size_um === missing
    @test any(d -> d.field == "pixel_size_um", disc4)
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: FAIL (`UndefVarError: derive_geometry not defined`).

- [ ] **Step 3: Append `derive_geometry` to `geometry.jl`**

```julia
# ---------------------------------------------------------------------------
# Geometry derivation orchestrator
# ---------------------------------------------------------------------------

"""
    GeometryDiscrepancy

A detected variance in a constant PRP field across multiple exposures,
or an unresolvable field (unknown detector). Surfaced as a Configuration-tab
banner (spec §6).
"""
struct GeometryDiscrepancy
    field   ::String
    message ::String
end

"""
    derive_geometry(prp_paths, setup_files) -> (geometry, discrepancies)

Sample the PRP files at `prp_paths` (reads all of them; at SSRL scale this is
~700 files × ~50 lines each, sub-second), derive per-field geometry with source
tags, and collect discrepancies.

`setup_files` is a vector of `analysis/setup_info_*.txt` paths; the latest by
filename sort is used for beam center + calibrated distance (filenames encode a
`YYYYMMDD_HHMMSS` timestamp so lexicographic sort is chronological).

Returns:
  `geometry` — a NamedTuple with fields:
      energy_kev, energy_kev_source,
      flight_path_m, flight_path_m_source,
      beam_center_x, beam_center_x_source,
      beam_center_y, beam_center_y_source,
      pixel_size_um, pixel_size_um_source
  `discrepancies` — a Vector{GeometryDiscrepancy}
"""
function derive_geometry(
    prp_paths::AbstractVector{<:AbstractString},
    setup_files::AbstractVector{<:AbstractString},
)
    discrepancies = GeometryDiscrepancy[]

    # Parse all PRPs (small files, sequential read is fine)
    parsed = [parse_prp(p) for p in prp_paths]

    # --- Energy (from PRP, should be constant) ---
    energy_vals = filter(!ismissing, [p.energy_kev for p in parsed])
    energy_kev = isempty(energy_vals) ? missing : first(energy_vals)
    if length(unique(round.(skipmissing([p.beam_energy_ev for p in parsed]); digits=2))) > 1
        push!(discrepancies, GeometryDiscrepancy("beam_energy_ev",
            "beam energy varies across PRPs: " *
            join(unique(round.(skipmissing([p.beam_energy_ev for p in parsed]); digits=2)), ", ") * " eV"))
    end

    # --- Pixel size from detector model (constant across run) ---
    detectors = unique(filter(!ismissing, [p.detector for p in parsed]))
    pixel_size_um  = missing
    pixel_size_source = "default"
    if length(detectors) == 1
        pitch = detector_pixel_size_um(detectors[1])
        if pitch === missing
            push!(discrepancies, GeometryDiscrepancy("pixel_size_um",
                "unknown detector model '$(detectors[1])': pixel pitch unknown, manual entry required"))
        else
            pixel_size_um    = pitch
            pixel_size_source = "prp"
        end
    elseif length(detectors) > 1
        push!(discrepancies, GeometryDiscrepancy("pixel_size_um",
            "multiple detector models found: $(join(detectors, ", "))"))
    end

    # --- Beam center + flight path from setup file ---
    beam_center_x = missing; beam_center_x_source = "default"
    beam_center_y = missing; beam_center_y_source = "default"
    flight_path_m = missing; flight_path_m_source = "default"

    if !isempty(setup_files)
        # Filenames are "setup_info_YYYYMMDD_HHMMSS.txt"; lexicographic sort picks the latest.
        latest_setup = last(sort(setup_files))
        si = parse_setup_info(latest_setup)
        if !ismissing(si.beam_center_x)
            beam_center_x = si.beam_center_x; beam_center_x_source = "setup"
        end
        if !ismissing(si.beam_center_y)
            beam_center_y = si.beam_center_y; beam_center_y_source = "setup"
        end
        if !ismissing(si.mean_distance_m)
            flight_path_m = si.mean_distance_m; flight_path_m_source = "setup"
        end
    end

    # Fallback: use PRP Pipe length when setup file absent or mean distance unparseable
    if ismissing(flight_path_m)
        pipe_vals = filter(!ismissing, [p.pipe_length_m for p in parsed])
        if !isempty(pipe_vals)
            flight_path_m = first(pipe_vals); flight_path_m_source = "prp"
            # Flag if pipe lengths vary (multi-setup discrepancy)
            if length(unique(round.(pipe_vals; digits=4))) > 1
                push!(discrepancies, GeometryDiscrepancy("flight_path_m",
                    "PRP pipe lengths vary: " * join(unique(pipe_vals .* 1000), ", ") * " mm"))
            end
        end
    end

    # Discrepancy: large gap between PRP pipe length and setup calibrated distance
    if !ismissing(flight_path_m) && flight_path_m_source == "setup"
        pipe_vals = filter(!ismissing, [p.pipe_length_m for p in parsed])
        if !isempty(pipe_vals)
            nominal = first(pipe_vals)
            frac_gap = abs(flight_path_m - nominal) / nominal
            if frac_gap > 0.01  # >1% gap (6.4% in the real data = 1809.5 vs 1700 mm)
                push!(discrepancies, GeometryDiscrepancy("flight_path_m_nominal_vs_calibrated",
                    "PRP pipe length $(round(nominal*1000; digits=1)) mm vs " *
                    "setup calibrated $(round(flight_path_m*1000; digits=1)) mm " *
                    "($(round(frac_gap*100; digits=1))% gap); using calibrated value"))
            end
        end
    end

    geometry = (
        energy_kev             = energy_kev,
        energy_kev_source      = ismissing(energy_kev) ? "default" : "prp",
        flight_path_m          = flight_path_m,
        flight_path_m_source   = flight_path_m_source,
        beam_center_x          = beam_center_x,
        beam_center_x_source   = beam_center_x_source,
        beam_center_y          = beam_center_y,
        beam_center_y_source   = beam_center_y_source,
        pixel_size_um          = pixel_size_um,
        pixel_size_um_source   = pixel_size_source,
    )
    return (geometry, discrepancies)
end
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/geometry.jl packages/HimalayaUI/test/test_ingestion_core.jl
git commit -m "feat(geometry): derive_geometry orchestrator with source tags and discrepancy detection"
```

---

## Task 4: `grouping.jl` — `scan_directory`

Enumerate TIF + PRP pairs from a directory using the existing `resolve_files` / `resolve_file_path` machinery from `config.jl`. Per spec §9.1 the dispatch argument to `resolve_files` is `ExperimentConfig`; this task addresses the "Directory-scan dispatch gap" open question by constructing a minimal config that is sufficient for enumeration.

**Files:**
- Create: `packages/HimalayaUI/src/grouping.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_core.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "scan_directory" begin
    dir = mktempdir()
    data_dir = joinpath(dir, "data")
    analysis_dir = joinpath(dir, "analysis")
    mkpath(data_dir); mkpath(analysis_dir)

    # Write 6 PRP files with 3 distinct horizontal positions (2 frames per slot)
    stems = ["HA_01_S001_0_001", "HA_01_S002_0_001",   # slot 1 (H ≈ 58.9)
             "HA_02_S003_0_001", "HA_02_S004_0_001",   # slot 2 (H ≈ 63.1)
             "HA_03_S005_0_001", "HA_03_S006_0_001"]   # slot 3 (H ≈ 67.3)
    h_positions = [58.9, 58.91, 63.1, 63.09, 67.3, 67.31]
    for (stem, hpos) in zip(stems, h_positions)
        write_prp(joinpath(data_dir, "$stem.prp");
            horizontal_position_mm = hpos)
        write(joinpath(data_dir, "$stem.tif"), "fake tif")
        write(joinpath(analysis_dir, "$stem.dat"), "fake dat")
    end

    # Write one setup_info file
    write_setup_info(joinpath(analysis_dir, "setup_info_20260425_181705.txt"))

    metas = HimalayaUI.scan_directory(data_dir, analysis_dir)

    # All 6 stems found
    @test length(metas) == 6
    # Each has a prp_path, tif_path, dat_path, and parsed prp fields
    m = first(metas)
    @test haskey(m, :stem)
    @test haskey(m, :prp_path)
    @test haskey(m, :tif_path)
    @test haskey(m, :prp)  # the parse_prp result
    # Stems returned sorted
    @test [m.stem for m in metas] == sort([m.stem for m in metas])
    # A missing PRP (has tif but no prp) → prp is missing
    write(joinpath(data_dir, "HA_04_S007_0_001.tif"), "fake tif")
    metas2 = HimalayaUI.scan_directory(data_dir, analysis_dir)
    @test length(metas2) == 7
    tif_only = first(filter(m -> m.stem == "HA_04_S007_0_001", metas2))
    @test tif_only.prp_path === nothing
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: FAIL.

- [ ] **Step 3: Create `grouping.jl` with `scan_directory`**

```julia
# packages/HimalayaUI/src/grouping.jl
#
# Directory scanning + sample grouping (spec §5, §9.1).
#
# scan_directory: enumerate (tif, prp, dat) triplets from a beamtime directory.
# group_into_samples: cluster exposures into loads and samples (Task 5–7).
# scan_and_group!: transactional orchestrator (ingest.jl, Task 8).
#
# Naming note: `auto_group` already exists in pipeline.jl for index-candidate
# grouping (an unrelated concept). This module uses `group_into_samples` (spec §9.1).

using Statistics

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
sidecar. Returns one `ExposureMeta` per TIFF stem, sorted by stem.

Reuses the `resolve_files` / `resolve_file_path` logic from `config.jl` for
consistent glob semantics, constructing a minimal `ExperimentConfig` for
dispatch (resolves the §9.1 open question: loosen dispatch vs. construct minimal
config — we construct a minimal config since `resolve_files` takes `ExperimentConfig`
for dispatch only and none of its fields are read by the dispatch branch).
"""
function scan_directory(
    data_dir::AbstractString,
    analysis_dir::AbstractString;
    tif_pattern::String = "{name}.tif",
    prp_pattern::String = "{name}.prp",
    dat_pattern::String = "{name}.dat",
)::Vector{ExposureMeta}
    # Build a minimal ExperimentConfig for the dispatch arg to resolve_files.
    # The dispatch method only needs the struct type; no fields are read.
    cfg = _minimal_scan_config(tif_pattern, prp_pattern, dat_pattern, data_dir, analysis_dir)

    # Enumerate all TIFs in data_dir using the empty prefix ""
    # (resolve_files("", ...) returns every file matching the pattern suffix).
    tif_stems = resolve_files(cfg, data_dir, "", tif_pattern)

    metas = ExposureMeta[]
    for stem in tif_stems
        tif_path = resolve_file_path(cfg, data_dir, stem, tif_pattern)
        prp_path = resolve_file_path(cfg, data_dir, stem, prp_pattern)
        dat_path = resolve_file_path(cfg, analysis_dir, stem, dat_pattern)
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
"""
function _minimal_scan_config(
    tif_pattern::String,
    prp_pattern::String,
    dat_pattern::String,
    data_dir::String,
    analysis_dir::String,
)::ExperimentConfig
    # ExperimentConfig field order confirmed from config.jl struct definition.
    # All manifest column fields default to 1 (unused); beamline fields to nothing.
    ExperimentConfig(
        # [experiment]
        "", "", "",
        # [beamline]
        nothing, nothing, "A-1", nothing, nothing, nothing,
        # [manifest]
        ",", 0, 0, 1, 1, 1, 1, 1, 1,
        # [layout]
        data_dir, analysis_dir, "simple",
        # [files]
        dat_pattern, tif_pattern,
    )
end
```

- [ ] **Step 4: Add to HimalayaUI.jl**

Insert after `include("geometry.jl")`:

```julia
include("geometry.jl")
include("grouping.jl")   # ← new
include("config.jl")
```

> **Note (verified against live `config.jl` 2026-06-18):** `ExperimentConfig` is a **plain `struct`** (`config.jl:21`) with **no `@kwdef` and no keyword constructor** — only the auto-generated positional constructor exists, so a keyword-form call (`ExperimentConfig(name=…, …)`) would `MethodError`. The 23-field positional order is, exactly: `name, description, manifest_file` (`[experiment]`); `energy_kev, flight_path_m, q_units, beam_center_x, beam_center_y, pixel_size_um` (`[beamline]`); `delimiter, skip_rows, header_row, col_sample_id, col_name, col_display_name, col_filenames, col_notes_sample, col_notes_exposure` (`[manifest]`); `data_dir, analysis_dir, exposure_type` (`[layout]`); `integration_pattern, image_pattern` (`[files]`). The `_minimal_scan_config` call above already matches this order exactly (the section comments line up 3/6/9/3/2). If a future field is added to the struct, this positional call must be updated in lockstep — re-read `config.jl:21-50` before editing.

- [ ] **Step 5: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/grouping.jl packages/HimalayaUI/src/HimalayaUI.jl packages/HimalayaUI/test/test_ingestion_core.jl
git commit -m "feat(grouping): scan_directory enumerates tif+prp+dat triplets"
```

---

## Task 5: `grouping.jl` — time-gap load segmentation

Implement `_segment_loads(metas)` — the first backbone step of §5: split the sorted-by-timestamp exposures into **loads** using time-gap detection. Gap-relative per spec §5: `gap > k × median(all_gaps)` with `k` defaulting to `10.0` (a clean bimodal histogram; 19 s intra-load vs 2578 s reload in the real data gives a ratio > 100). Includes the unimodal fallback (no clear bimodal separation → one load, flag).

**Files:**
- Modify: `packages/HimalayaUI/src/grouping.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_core.jl`

- [ ] **Step 1: Write the failing test**

The test uses synthetic timestamps that mirror the real data pattern (19 s intra-load, 2578 s reload gap between two loads of 4 exposures each):

```julia
@testset "load segmentation" begin
    using Dates

    # Build 8 ExposureMeta with synthetic timestamps:
    #   4 in Load 1 at T+0, T+19, T+38, T+57
    #   2578s gap (reload)
    #   4 in Load 2 at T+2635, T+2654, T+2673, T+2692
    t0 = DateTime(2026, 4, 26, 23, 14, 8)
    offsets = [0, 19, 38, 57, 2635, 2654, 2673, 2692]  # seconds
    metas = [ExposureMeta(
        "HA_$(lpad(i,2,'0'))_S$(lpad(i,4,'0'))_0_001",
        nothing, nothing, nothing,
        (timestamp = t0 + Second(s), horizontal_position_mm = 58.9 + Float64(i)*4.0,
         beam_energy_ev = 9000.0, energy_kev = 9.0, pipe_length_m = 1.7,
         detector = "Pilatus 1M", exposure_time_s = 15.0),
    ) for (i, s) in enumerate(offsets)]

    loads = HimalayaUI._segment_loads(metas)
    @test length(loads) == 2
    @test length(loads[1]) == 4
    @test length(loads[2]) == 4
    @test loads[1][1].stem == "HA_01_S0001_0_001"
    @test loads[2][1].stem == "HA_05_S0005_0_001"

    # Unimodal fallback: all gaps equal → one load, returned with flag
    uniform_metas = [ExposureMeta(
        "HA_$(lpad(i,2,'0'))_S$(lpad(i,4,'0'))_0_001",
        nothing, nothing, nothing,
        (timestamp = t0 + Second(i * 20), horizontal_position_mm = 58.9 + Float64(i)*4.0,
         beam_energy_ev = 9000.0, energy_kev = 9.0, pipe_length_m = 1.7,
         detector = "Pilatus 1M", exposure_time_s = 15.0),
    ) for i in 1:4]
    uni_loads, uni_flag = HimalayaUI._segment_loads_with_flag(uniform_metas)
    @test length(uni_loads) == 1
    @test uni_flag == :unimodal_fallback

    # Single exposure → one load, no-crash
    single = [first(metas)]
    sl, sf = HimalayaUI._segment_loads_with_flag(single)
    @test length(sl) == 1
    @test sf == :ok
end
```

> **Note on ExposureMeta constructor:** the `prp` field is typed `Union{NamedTuple, Nothing}`. The test passes a NamedTuple literal directly; this works because NamedTuple is a concrete type. If the struct requires a specific named NamedTuple type, wrap in a named-tuple literal: `(timestamp=…, horizontal_position_mm=…, …)`.

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: FAIL.

- [ ] **Step 3: Append load segmentation to `grouping.jl`**

```julia
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

`flag ∈ {:ok, :unimodal_fallback, :single_exposure}`:
- `:ok`                — at least one clear bimodal split found
- `:unimodal_fallback` — no gap exceeds the threshold; entire directory is one load
- `:single_exposure`   — only one exposure; trivially one load
"""
function _segment_loads_with_flag(metas::Vector{ExposureMeta}; gap_k::Float64 = 10.0)
    isempty(metas) && return (Vector{ExposureMeta}[], :ok)

    # Sort by timestamp; push exposures with missing timestamp to the end.
    sorted = sort(metas; by = m -> begin
        ts = m.prp !== nothing ? m.prp.timestamp : missing
        ismissing(ts) ? DateTime(9999) : ts
    end)

    length(sorted) == 1 && return ([[sorted[1]]], :single_exposure)

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

    med = Statistics.median(gaps)
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
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/grouping.jl packages/HimalayaUI/test/test_ingestion_core.jl
git commit -m "feat(grouping): time-gap load segmentation with unimodal fallback"
```

---

## Task 6: `grouping.jl` — stepping-axis detection + slot clustering

Implement `_detect_stepping_axis(load_metas)` and `_cluster_slots(load_metas)` — spec §5 step 1: the motor with **dominant inter-burst variance** is the stepping axis; position-gap clustering groups frames into slots. The real data uses `horizontal_position_mm`; the code detects this rather than hard-coding it.

**Files:**
- Modify: `packages/HimalayaUI/src/grouping.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_core.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "slot clustering" begin
    using Dates

    # Simulate one load: 3 slots × 4 frames each, H positions stepped ~4 mm apart,
    # within-burst jitter ≈ 0.3 mm. Mirrors the real HA data pattern from spec §5.
    t0 = DateTime(2026, 4, 26, 23, 14, 8)
    slot_h = [74.80, 70.85, 67.22]  # three slot centers
    jitter = [0.0, 0.13, -0.24, 0.05]  # within-burst jitter, 4 frames
    frames = ExposureMeta[]
    for (si, h_center) in enumerate(slot_h)
        for (fi, j) in enumerate(jitter)
            offset = (si - 1) * 4 * 19 + (fi - 1) * 19
            stem = "HA_$(si)_$(fi)_S$(lpad((si-1)*4+fi, 4,'0'))_0_001"
            push!(frames, ExposureMeta(stem, nothing, nothing, nothing,
                (timestamp = t0 + Second(offset),
                 horizontal_position_mm = h_center + j,
                 beam_energy_ev = 9000.0, energy_kev = 9.0, pipe_length_m = 1.7,
                 detector = "Pilatus 1M", exposure_time_s = 15.0)))
        end
    end

    slots = HimalayaUI._cluster_slots(frames)
    @test length(slots) == 3
    @test all(length(s) == 4 for s in slots)

    # All 4 frames of slot 1 should have H ≈ 74.80
    slot1_h = [m.prp.horizontal_position_mm for m in slots[1]]
    @test all(abs(h - 74.80) < 0.5 for h in slot1_h)

    # Median-delta-near-zero fallback (the `med_delta < 1e-6` branch).
    #
    # This exercises the documented algorithm intent: within-slot revisits collapse the
    # MEDIAN consecutive-frame delta to ~0 (most consecutive frames sit at the *same*
    # position), so the plain `med_delta × slot_k` tolerance is degenerate (0) and would
    # split on every frame. The fallback instead learns the slot-spacing tolerance from
    # the non-zero deltas (the burst→burst jumps). Fixture: 3 slots, each a 3-frame burst
    # at a fixed position (consecutive in-burst deltas == 0), with ~4 mm jumps between
    # bursts. The 6 in-burst deltas are 0 and the 2 between-burst deltas are ~4 mm, so
    # median(all 8) == 0 → fallback branch → tolerance from median(nonzero)=~4 mm × 5.
    burst_frames = ExposureMeta[]
    for (si, h_center) in enumerate([74.80, 70.85, 67.22])
        for _ in 1:3
            push!(burst_frames, ExposureMeta("b_$(si)_$(length(burst_frames))",
                nothing, nothing, nothing,
                (timestamp = t0 + Second(length(burst_frames) * 19),
                 horizontal_position_mm = h_center,   # identical within the burst → zero deltas
                 beam_energy_ev = 9000.0, energy_kev = 9.0, pipe_length_m = 1.7,
                 detector = "Pilatus 1M", exposure_time_s = 15.0)))
        end
    end

    slots_bf = HimalayaUI._cluster_slots(burst_frames)
    @test length(slots_bf) == 3                       # one slot per burst position
    @test all(length(s) == 3 for s in slots_bf)       # 3 frames each
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: FAIL.

- [ ] **Step 3: Append slot clustering to `grouping.jl`**

```julia
# ---------------------------------------------------------------------------
# Stepping-axis detection + slot clustering (spec §5, step 1)
# ---------------------------------------------------------------------------

"""
    _cluster_slots(load_metas; slot_k = 5.0) -> Vector{Vector{ExposureMeta}}

Group the exposures in one load into slots (one slot = one sample position).

Algorithm (spec §5):
1. Extract horizontal positions; if all are missing → one slot (fallback).
2. Compute consecutive-pair position deltas; the within-burst jitter tolerance
   is `slot_k × median(|consecutive deltas|)`.
3. A position gap exceeding that tolerance starts a new slot.
4. Median-near-zero fallback (spec §5): when the *median* consecutive delta is ~0
   — i.e. most consecutive frames sit at the *same* position (multi-frame bursts /
   within-slot revisits) — the plain `median × slot_k` tolerance degenerates to 0
   and would split on every step. Instead learn the tolerance from the non-zero
   deltas (the burst→burst jumps): `median(nonzero deltas) × slot_k`. If there are
   no non-zero deltas at all (every frame at one position), tolerance = `Inf` → one
   slot.

`slot_k = 5.0` is a validated default (real data: 0.30 mm within-burst jitter
vs 3.95 mm slot spacing → ratio 13; a 5× multiplier sits well inside the gap).

KNOWN LIMITATION (deferred — see spec §5 "single-frame acquisitions … else flag"):
a load of pure *single-frame* acquisitions visited round-robin (each slot exactly
one frame, no multi-frame burst anywhere) has every consecutive delta equal to a
slot-to-slot jump, so the jitter population and the spacing population are
indistinguishable — there is nothing to "learn the gap from". This case is NOT
handled here (it would under-split to one slot); the spec's prescribed behavior is
to flag for human review. Wiring that flag is deferred to Phase D's grouping-review
UI; for Phase B it falls through to a single slot. See the Self-Review note.
"""
function _cluster_slots(
    load_metas::Vector{ExposureMeta};
    slot_k::Float64 = 5.0,
)::Vector{Vector{ExposureMeta}}
    isempty(load_metas) && return Vector{ExposureMeta}[]
    length(load_metas) == 1 && return [[load_metas[1]]]

    hpos = [m.prp !== nothing && !ismissing(m.prp.horizontal_position_mm) ?
             m.prp.horizontal_position_mm : missing
            for m in load_metas]

    all_missing = all(ismissing, hpos)
    if all_missing
        # No position axis available: one slot, flag for review
        return [load_metas]
    end

    # Consecutive position deltas
    deltas = Float64[]
    for i in 2:length(hpos)
        if !ismissing(hpos[i]) && !ismissing(hpos[i-1])
            push!(deltas, abs(hpos[i] - hpos[i-1]))
        end
    end

    if isempty(deltas)
        return [load_metas]
    end

    med_delta = Statistics.median(deltas)
    # Median-near-zero fallback: when most consecutive frames sit at the same position
    # (multi-frame bursts), the median delta is ~0 and `median × slot_k` would be a
    # degenerate (zero) tolerance. Learn the slot-spacing tolerance from the non-zero
    # deltas (the burst→burst jumps) instead; Inf when there are no jumps at all.
    # NOTE: in the fallback regime the non-zero deltas ARE the inter-slot spacing, so the
    # tolerance must sit BELOW them — divide by slot_k, not multiply (spec §5: tolerance
    # between jitter and slot spacing). The normal branch multiplies because there
    # median(deltas) is the within-burst JITTER. (Corrected 2026-06-19, Phase-B Task 6.)
    local threshold::Float64
    if med_delta < 1e-6
        nonzero = filter(d -> d > 1e-6, deltas)
        threshold = isempty(nonzero) ? Inf : Statistics.median(nonzero) / slot_k
    else
        threshold = med_delta * slot_k
    end

    slots = Vector{ExposureMeta}[[load_metas[1]]]
    for i in 2:length(load_metas)
        h_prev = ismissing(hpos[i-1]) ? missing : hpos[i-1]
        h_curr = ismissing(hpos[i])   ? missing : hpos[i]
        gap = (!ismissing(h_prev) && !ismissing(h_curr)) ? abs(h_curr - h_prev) : 0.0
        if gap > threshold
            push!(slots, ExposureMeta[])
        end
        push!(last(slots), load_metas[i])
    end
    return slots
end
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/grouping.jl packages/HimalayaUI/test/test_ingestion_core.jl
git commit -m "feat(grouping): stepping-axis detection and slot clustering"
```

---

## Task 7: `grouping.jl` — `group_into_samples` + auto-naming

Assemble the load/slot structure into concrete `GroupedLoad` / `GroupedSample` / `GroupedExposure` structs that `scan_and_group!` will persist. Includes the auto-naming rule from spec §5: `<label> (SNNPMM)` where S=load index and P=slot index within load.

**Files:**
- Modify: `packages/HimalayaUI/src/grouping.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_core.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "group_into_samples" begin
    using Dates

    # Two-load synthetic directory: 2 loads × 3 slots × 2 frames
    t0 = DateTime(2026, 4, 26, 23, 14, 8)
    slot_h_l1 = [74.80, 70.85, 67.22]
    slot_h_l2 = [74.75, 70.90, 67.18]  # same positions in load 2 (reload)
    frames = ExposureMeta[]
    global_s = 1
    for (li, slot_h) in enumerate([slot_h_l1, slot_h_l2])
        load_start = li == 1 ? 0 : 600  # 600 s reload gap
        for (si, h) in enumerate(slot_h)
            for fi in 1:2
                offset = load_start + (si - 1) * 2 * 19 + (fi - 1) * 19
                stem = "HA_$(li)_$(si)_$(fi)_S$(lpad(global_s, 4,'0'))_0_001"
                global_s += 1
                push!(frames, ExposureMeta(stem, "/d/$stem.tif", "/d/$stem.prp", "/a/$stem.dat",
                    (timestamp = t0 + Second(offset),
                     horizontal_position_mm = h + (fi == 2 ? 0.1 : 0.0),
                     beam_energy_ev = 9000.0, energy_kev = 9.0, pipe_length_m = 1.7,
                     detector = "Pilatus 1M", exposure_time_s = 15.0)))
            end
        end
    end

    result = HimalayaUI.group_into_samples(frames)

    # Structure: 2 loads
    @test length(result.loads) == 2
    # 3 slots per load
    @test all(length(l.samples) == 3 for l in result.loads)
    # 2 frames per sample
    @test all(length(s.exposures) == 2 for l in result.loads for s in l.samples)

    # Auto-naming: first sample of load 1, slot 1 → "HA ... (S01P01)"
    name11 = result.loads[1].samples[1].name
    @test occursin("S01P01", name11)
    # Second load: S02P01
    name21 = result.loads[2].samples[1].name
    @test occursin("S02P01", name21)

    # Load indices 1-based
    @test result.loads[1].load_index == 1
    @test result.loads[2].load_index == 2

    # Frame counts
    @test result.loads[1].frame_count == 6  # 3 slots × 2 frames
    @test result.loads[2].frame_count == 6

    # Load time range
    @test result.loads[1].start_time isa DateTime
    @test result.loads[1].end_time   isa DateTime
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: FAIL.

- [ ] **Step 3: Append the grouped structs + `group_into_samples` to `grouping.jl`**

```julia
# ---------------------------------------------------------------------------
# Grouped structs returned by group_into_samples
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
P=slot index (both 1-based, zero-padded). `label_hint` is derived from the
filename token when parseable; empty string yields just the coordinate anchor.
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
  1. Sort by timestamp.
  2. Time-gap segment into loads (_segment_loads_with_flag).
  3. Within each load, cluster by position (_cluster_slots).
  4. Build GroupedLoad/GroupedSample/GroupedExposure tree with auto-names.
"""
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
            # Use the most common label hint across the slot's stems
            label_hint = _label_from_stem(first(slot_metas).stem)

            exps = [begin
                sid, fno = _parse_scan_frame(m.stem)
                GroupedExposure(
                    m.stem,
                    m.tif_path, m.prp_path, m.dat_path,
                    m.prp !== nothing && !ismissing(m.prp.timestamp) ? m.prp.timestamp : missing,
                    m.prp !== nothing ? m.prp.exposure_time_s : missing,
                    m.prp !== nothing ? m.prp.horizontal_position_mm : missing,
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

        timestamps = filter(!ismissing, [begin
            ts = m.prp !== nothing ? m.prp.timestamp : missing
            ismissing(ts) ? missing : ts
        end for m in load_metas])
        start_time = isempty(timestamps) ? missing : minimum(timestamps)
        end_time   = isempty(timestamps) ? missing : maximum(timestamps)

        push!(grouped_loads, GroupedLoad(
            li,
            length(load_metas),
            start_time,
            end_time,
            grouped_samples,
            seg_flag,   # :ok | :unimodal_fallback | :single_exposure — propagated from segmentation (spec §5)
        ))
    end

    return GroupingResult(grouped_loads, discrepancies)
end
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/grouping.jl packages/HimalayaUI/test/test_ingestion_core.jl
git commit -m "feat(grouping): group_into_samples produces GroupedLoad/Sample/Exposure tree"
```

---

## Task 8: `ingest.jl` — `scan_and_group!` transactional orchestrator

The DB-writing entry point. Creates (or updates) the experiment's geometry, creates `loads`, `samples`, and `exposures` rows via the Phase-A writers, then fires `analyze_exposure!` for each new exposure. Clones `_reingest_inner!`'s insert-only discipline and wraps the whole scan in `_DB_WRITE_LOCK`. The Phase-A `create_exposure!` already accepts `sample_id = nothing` (transient scan state), so analysis can run after grouping is complete.

**Files:**
- Create: `packages/HimalayaUI/src/ingest.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl`
- Test: `packages/HimalayaUI/test/test_ingestion_core.jl`

- [ ] **Step 1: Write the failing test**

This test exercises the full `scan_and_group!` path without real `.dat` files. It stubs analysis by verifying that exposures land in the DB with the correct `experiment_id` and that `sample_id` is populated after grouping.

```julia
@testset "scan_and_group! inserts loads/samples/exposures" begin
    dir = mktempdir()
    data_dir     = joinpath(dir, "data")
    analysis_dir = joinpath(dir, "analysis")
    mkpath(data_dir); mkpath(analysis_dir)

    # 2 slots × 2 frames in one load
    slot_h = [58.9, 63.1]
    t0 = DateTime(2026, 4, 26, 23, 14, 8)
    all_stems = String[]
    for (si, h) in enumerate(slot_h)
        for fi in 1:2
            offset = (si - 1) * 2 * 19 + (fi - 1) * 19
            stem = "HA_$(si)_$(fi)_S$(lpad((si-1)*2+fi, 4,'0'))_0_001"
            push!(all_stems, stem)
            write_prp(joinpath(data_dir, "$stem.prp");
                timestamp = Dates.format(t0 + Second(offset), "dd u yyyy HH:MM:SS"),
                beam_energy_ev = 9000.027604502573,
                pipe_length_mm = 1700, detector = "Pilatus 1M",
                exposure_time = 15.0, horizontal_position_mm = h)
            write(joinpath(data_dir, "$stem.tif"), "fake tif")
            # No .dat files: analyze_exposure! will fail to read them, which is OK
            # for this insert-path test (we set analyze=false)
        end
    end
    write_setup_info(joinpath(analysis_dir, "setup_info_20260425_181705.txt"))

    db = fresh_db()
    # Create a bare experiment (Phase-A create_experiment! signature)
    exp_id = HimalayaUI.create_experiment!(db;
        name = "test-ingest", path = dir,
        data_dir = data_dir, analysis_dir = analysis_dir)

    HimalayaUI.scan_and_group!(db, exp_id, dir; analyze = false)

    # Exposures in DB
    exps = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, experiment_id, sample_id, filename, image_path, scan_id, frame_no FROM exposures WHERE experiment_id = ?",
        [exp_id]))
    @test length(exps) == 4
    @test all(e.experiment_id == exp_id for e in exps)
    # After grouping, sample_id is populated
    @test all(!ismissing(e.sample_id) for e in exps)
    # image_path persisted (non-NULL, equals the written .tif path) — the image route
    # and prewarm_thumbnails! both filter `WHERE image_path IS NOT NULL`.
    @test all(!ismissing(e.image_path) for e in exps)
    @test all(endswith(String(e.image_path), ".tif") for e in exps)
    # scan_id / frame_no parsed from the filename stem (`_S<digits>_…_<frame>`) and persisted (spec §4/§10)
    @test all(!ismissing(e.scan_id) for e in exps)
    @test all(e.frame_no == 1 for e in exps)   # every fixture stem ends `_001`

    # Samples: 2 (one per slot)
    samps = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, name, load_id, slot_index FROM samples WHERE experiment_id = ?", [exp_id]))
    @test length(samps) == 2
    @test all(s.slot_index in [1, 2] for s in samps)
    @test any(occursin("S01P01", s.name) for s in samps)

    # Loads: 1
    loads = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, load_index, frame_count FROM loads WHERE experiment_id = ?", [exp_id]))
    @test length(loads) == 1
    @test loads[1].load_index == 1
    @test loads[1].frame_count == 4

    # Geometry written to experiments row
    row = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT energy_kev, flight_path_m, flight_path_m_source, beam_center_x FROM experiments WHERE id = ?",
        [exp_id])))
    @test row.energy_kev ≈ 9.000027604502573
    @test row.flight_path_m ≈ 1.8095   # from setup file
    @test row.flight_path_m_source == "setup"
    @test row.beam_center_x ≈ 421.409

    # Idempotent: a second scan of the same directory is a true no-op. Insert-only
    # discipline applies to loads AND samples AND exposures — re-running must not mint
    # a fresh load_id (which would re-key the sample dedup and duplicate samples).
    HimalayaUI.scan_and_group!(db, exp_id, dir; analyze = false)
    exps2 = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM exposures WHERE experiment_id = ?", [exp_id]))
    @test length(exps2) == 4  # exposures unchanged
    loads2 = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM loads WHERE experiment_id = ?", [exp_id]))
    @test length(loads2) == 1  # loads unchanged (no orphan duplicate load row)
    samps2 = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM samples WHERE experiment_id = ?", [exp_id]))
    @test length(samps2) == 2  # samples unchanged (sample dedup keyed on a stable load_id)
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: FAIL (`UndefVarError: scan_and_group! not defined`).

- [ ] **Step 3: Create `ingest.jl`**

```julia
# packages/HimalayaUI/src/ingest.jl
#
# scan_and_group!: the transactional ingest orchestrator (spec §9.1).
#
# Discipline cloned from _reingest_inner! (cli.jl:156):
#   - Always insert-only for exposures (dedup key (experiment_id, filename)).
#   - Never clobber human-set fields (name_source='human', spec §4).
#   - The DB write is wrapped in _DB_WRITE_LOCK (server.jl:29).
#   - analyze_exposure! is called OUTSIDE the main transaction (like cli_init_with_db!).
#
# Note: `_reingest_inner!` is not called through; scan_and_group! is net-new and
# clones the discipline because _reingest_inner! requires a manifest (returns
# :no_manifest early) and is insert-only for exposures but clobbers sample fields.

"""
    scan_and_group!(db, experiment_id, root_dir; additive=true, analyze=true,
                    tif_pattern="{name}.tif", prp_pattern="{name}.prp", dat_pattern="{name}.dat")

Full ingest of a beamtime directory into `db` under `experiment_id`.

Steps:
  1. Resolve `data_dir` and `analysis_dir` from the `experiments` row.
  2. Scan: enumerate TIF+PRP+DAT triplets with `scan_directory`.
  3. Geometry: derive from the sampled PRPs + the latest setup_info file; write
     to `experiments` row (respects `*_source='human'` never-clobber).
  4. Group: `group_into_samples` → loads/samples/slots.
  5. Persist: insert `loads`, `samples`, `exposures` rows (insert-only per the
     dedup key `(experiment_id, filename)` — no clobber on rescan).
  6. Analyze: run `analyze_exposure!` for every newly inserted exposure (skips
     those already in the DB). Runs OUTSIDE the write transaction (same contract
     as `cli_init_with_db!`).

When `additive=true` (the default for rescan), files already in the DB are
skipped silently (UNIQUE constraint on `(experiment_id, filename)` enforces
this at the DB level). When `analyze=false`, step 6 is skipped (useful for
tests and the HTTP "scan-only" path).
"""
function scan_and_group!(
    db         ::SQLite.DB,
    experiment_id::Int,
    root_dir   ::AbstractString;
    additive   ::Bool   = true,
    analyze    ::Bool   = true,
    tif_pattern::String = "{name}.tif",
    prp_pattern::String = "{name}.prp",
    dat_pattern::String = "{name}.dat",
)
    root_dir = abspath(root_dir)

    # Resolve data_dir and analysis_dir from the experiment row.
    exp_row = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT data_dir, analysis_dir FROM experiments WHERE id = ?", [experiment_id])))
    data_dir     = String(exp_row.data_dir)
    analysis_dir = String(exp_row.analysis_dir)

    # -----------------------------------------------------------------------
    # 1. Scan: enumerate TIF+PRP+DAT triplets
    # -----------------------------------------------------------------------
    metas = scan_directory(data_dir, analysis_dir;
        tif_pattern = tif_pattern,
        prp_pattern = prp_pattern,
        dat_pattern = dat_pattern)

    isempty(metas) && return (status = :empty, added_loads = 0, added_samples = 0, added_exposures = 0)

    # -----------------------------------------------------------------------
    # 2. Geometry: derive + write to experiments row (never-clobber human fields)
    # -----------------------------------------------------------------------
    prp_paths   = filter(!isnothing, [m.prp_path for m in metas])
    setup_files = filter(isfile, let
        isdir(analysis_dir) ?
            [joinpath(analysis_dir, f)
             for f in readdir(analysis_dir)
             if startswith(f, "setup_info_") && endswith(f, ".txt")] :
            String[]
    end)

    geo, _disc = derive_geometry(String.(prp_paths), String.(setup_files))

    # Persist geometry with never-clobber: only write fields whose current source
    # is not 'human' (spec §4).
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            _update_geometry_if_not_human!(db, experiment_id, geo)
        end
    end

    # -----------------------------------------------------------------------
    # 3. Group: loads / samples / slots
    # -----------------------------------------------------------------------
    grouping = group_into_samples(metas)

    # -----------------------------------------------------------------------
    # 4. Persist: loads → samples → exposures (all inside one write lock)
    # -----------------------------------------------------------------------
    new_exposure_ids = Int[]

    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            # Existing exposures (dedup key)
            existing_filenames = Set(String.(getproperty.(
                Tables.rowtable(DBInterface.execute(db,
                    "SELECT filename FROM exposures WHERE experiment_id = ?", [experiment_id])),
                :filename)))

            for gl in grouping.loads
                # Dedup the load row by (experiment_id, load_index). The loads table
                # has only a non-unique index (Phase A Task 2: loads_experiment_idx),
                # so an unconditional INSERT would mint a fresh load_id on every rescan
                # — that would re-key the sample dedup below onto a NEW load_id, find
                # nothing, and create duplicate samples + an orphan load row. Mirror the
                # sample/exposure insert-only discipline: reuse the existing load_id when
                # present, insert only when absent. This makes a clean rescan a true no-op.
                existing_load = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id FROM loads WHERE experiment_id = ? AND load_index = ?",
                    [experiment_id, gl.load_index]))
                if isempty(existing_load)
                    load_id = DBInterface.lastrowid(DBInterface.execute(db, """
                        INSERT INTO loads (experiment_id, load_index, frame_count,
                                           start_time, end_time)
                        VALUES (?, ?, ?, ?, ?)
                    """, [experiment_id, gl.load_index, gl.frame_count,
                          ismissing(gl.start_time) ? nothing : Dates.format(gl.start_time, "yyyy-mm-ddTHH:MM:SS"),
                          ismissing(gl.end_time)   ? nothing : Dates.format(gl.end_time,   "yyyy-mm-ddTHH:MM:SS")]))
                else
                    load_id = Int(first(existing_load).id)
                end

                for gs in gl.samples
                    # Create sample row (never-clobber: skip if (load_id, slot_index) exists)
                    existing_sample = Tables.rowtable(DBInterface.execute(db,
                        "SELECT id FROM samples WHERE experiment_id = ? AND load_id = ? AND slot_index = ?",
                        [experiment_id, load_id, gs.slot_index]))
                    if isempty(existing_sample)
                        sample_id = create_sample!(db;
                            experiment_id  = experiment_id,
                            name           = gs.name,
                            load_id        = Int(load_id),
                            slot_index     = gs.slot_index,
                            grouping_source = gs.grouping_source,
                            name_source    = gs.name_source)
                    else
                        sample_id = Int(first(existing_sample).id)
                    end

                    for ge in gs.exposures
                        ge.stem in existing_filenames && continue  # insert-only
                        exp_id = create_exposure!(db;
                            experiment_id          = experiment_id,
                            sample_id              = sample_id,
                            filename               = ge.stem,
                            image_path             = ge.tif_path,   # required by the image route + prewarm_thumbnails! (WHERE image_path IS NOT NULL)
                            prp_path               = ge.prp_path,
                            timestamp              = ismissing(ge.timestamp) ? nothing :
                                                     Dates.format(ge.timestamp, "yyyy-mm-ddTHH:MM:SS"),
                            exposure_time          = ismissing(ge.exposure_time_s) ? nothing :
                                                     ge.exposure_time_s,
                            horizontal_position    = ismissing(ge.horizontal_position_mm) ? nothing :
                                                     ge.horizontal_position_mm,
                            scan_id                = ismissing(ge.scan_id)  ? nothing : ge.scan_id,
                            frame_no               = ismissing(ge.frame_no) ? nothing : ge.frame_no,
                            load_id                = Int(load_id))
                        push!(new_exposure_ids, exp_id)
                        push!(existing_filenames, ge.stem)
                    end
                end
            end
        end
    end

    # -----------------------------------------------------------------------
    # 5. Analyze new exposures OUTSIDE the write transaction (same contract
    #    as cli_init_with_db!: a crash mid-analyze must not roll back ingest)
    # -----------------------------------------------------------------------
    if analyze
        for eid in new_exposure_ids
            try
                analyze_exposure!(db, eid)
            catch e
                @warn "scan_and_group!: analyze_exposure! failed" exposure_id=eid exception=e
            end
        end
    end

    return (
        status           = :ok,
        added_loads      = length(grouping.loads),
        added_samples    = sum(length(l.samples) for l in grouping.loads),
        added_exposures  = length(new_exposure_ids),
    )
end

"""
    _update_geometry_if_not_human!(db, experiment_id, geo)

Write derived geometry fields to the `experiments` row, skipping any field
whose current `*_source` is `'human'` (never-clobber contract, spec §4).
Called inside a write transaction.
"""
function _update_geometry_if_not_human!(db::SQLite.DB, experiment_id::Int, geo::NamedTuple)
    current = first(Tables.rowtable(DBInterface.execute(db,
        """SELECT energy_kev_source, flight_path_m_source,
                  beam_center_x_source, beam_center_y_source, pixel_size_um_source
             FROM experiments WHERE id = ?""", [experiment_id])))

    updates = Pair{String,Any}[]

    function maybe!(field, source_field, val, source)
        curr_src = getproperty(current, Symbol(source_field))
        if !ismissing(curr_src) && String(curr_src) == "human"
            return  # never overwrite a human-set field
        end
        val === missing && return
        push!(updates, field => val)
        push!(updates, source_field => source)
    end

    maybe!("energy_kev",    "energy_kev_source",    geo.energy_kev,    geo.energy_kev_source)
    maybe!("flight_path_m", "flight_path_m_source", geo.flight_path_m, geo.flight_path_m_source)
    maybe!("beam_center_x", "beam_center_x_source", geo.beam_center_x, geo.beam_center_x_source)
    maybe!("beam_center_y", "beam_center_y_source", geo.beam_center_y, geo.beam_center_y_source)
    maybe!("pixel_size_um", "pixel_size_um_source", geo.pixel_size_um, geo.pixel_size_um_source)

    isempty(updates) && return

    cols = join(first.(updates) .* " = ?", ", ")
    vals = last.(updates)
    DBInterface.execute(db,
        "UPDATE experiments SET $cols WHERE id = ?",
        vcat(vals, [experiment_id]))
end
```

- [ ] **Step 4: Add to HimalayaUI.jl**

Insert after `include("grouping.jl")` and before `include("config.jl")`:

```julia
include("grouping.jl")
include("ingest.jl")     # ← new
include("config.jl")
```

- [ ] **Step 5: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/ingest.jl packages/HimalayaUI/src/HimalayaUI.jl packages/HimalayaUI/test/test_ingestion_core.jl
git commit -m "feat(ingest): scan_and_group! transactional orchestrator"
```

---

## Task 9: Regression floor test against real SSRL fixtures

Add a test that runs `scan_directory` + `group_into_samples` against a **small synthetic fixture** that numerically mirrors the real SSRL 2026-04 run parameters (no dependency on `/Volumes/data`). Assert regression floors, not exact counts, per the repo convention ("regression floors, not hard-coded counts" from `CLAUDE.md`).

**Files:**
- Test: `packages/HimalayaUI/test/test_ingestion_core.jl`

- [ ] **Step 1: Write the test**

The fixture mirrors Load 1 of the 1p7m run: 4 slots of ~4 frames each, ~4 mm spacing, ~19 s within-burst gaps, preceded by a large pre-load gap from the calibration exposures (15867 s). The numeric parameters come from real data analysis conducted 2026-06-18.

```julia
@testset "grouping regression floor — SSRL 2026-04 Load 1 fixture" begin
    using Dates

    # Mirror the real 1p7m run structure (Load 1 of 13):
    # 4 slots, 4 frames each, H steps ~3.95 mm, within-burst jitter ≈ 0.3 mm,
    # 19 s between frames, 15 s exposure time (from PRP).
    # Confirmed parameters from real data sweep 2026-06-18.
    dir = mktempdir()
    data_dir     = joinpath(dir, "data")
    analysis_dir = joinpath(dir, "analysis")
    mkpath(data_dir); mkpath(analysis_dir)

    SLOT_CENTERS_MM = [74.80, 70.85, 67.22, 63.49]  # from real HA data
    WITHIN_BURST_JITTER = [0.0, 0.13, -0.24, 0.05]
    FRAME_GAP_S = 19
    t0 = DateTime(2026, 4, 25, 23, 14, 8)

    scan_id = 2001
    stems = String[]
    for (si, h_center) in enumerate(SLOT_CENTERS_MM)
        for (fi, jitter) in enumerate(WITHIN_BURST_JITTER)
            offset = (si - 1) * length(WITHIN_BURST_JITTER) * FRAME_GAP_S + (fi - 1) * FRAME_GAP_S
            stem = "HA_$(si)_$(lpad(scan_id, 4,'0'))_S$(scan_id)_0_001"
            scan_id += 1
            push!(stems, stem)
            write_prp(joinpath(data_dir, "$stem.prp");
                timestamp        = Dates.format(t0 + Second(offset), "dd u yyyy HH:MM:SS"),
                beam_energy_ev   = 9000.027604502573,
                pipe_length_mm   = 1700,
                detector         = "Pilatus 1M",
                exposure_time    = 15.0,
                horizontal_position_mm = h_center + jitter)
            write(joinpath(data_dir, "$stem.tif"), "fake tif")
        end
    end
    write_setup_info(joinpath(analysis_dir, "setup_info_20260425_181705.txt");
        beam_center_x = 421.409, beam_center_y = 836.946, mean_distance_mm = 1809.5)

    metas = HimalayaUI.scan_directory(data_dir, analysis_dir)

    # Regression floors for scan_directory
    @test length(metas) >= 16          # ≥ all 16 exposures found
    @test length(metas) <= 16          # ≤ only the 16 we wrote (no phantom)

    result = HimalayaUI.group_into_samples(metas)

    # One load (19 s gaps are well within one load)
    @test length(result.loads) >= 1
    # At least 4 slots detected
    total_samples = sum(length(l.samples) for l in result.loads)
    @test total_samples >= 4           # floor: found at least 4 slots
    @test total_samples <= 6           # ceiling: shouldn't massively over-split
    # Frame count matches
    total_frames = sum(l.frame_count for l in result.loads)
    @test total_frames == 16

    # Geometry derivation
    prp_paths   = filter(!isnothing, [m.prp_path for m in metas])
    setup_files = [joinpath(analysis_dir, "setup_info_20260425_181705.txt")]
    geo, disc   = HimalayaUI.derive_geometry(String.(prp_paths), String.(setup_files))

    @test geo.energy_kev ≈ 9.000027604502573
    @test geo.flight_path_m ≈ 1.8095            # calibrated, not PRP nominal
    @test geo.flight_path_m_source == "setup"
    @test geo.beam_center_x ≈ 421.409
    @test geo.pixel_size_um ≈ 172.0

    # The nominal vs calibrated gap is flagged as a discrepancy (6.4% gap)
    @test any(d -> occursin("flight_path_m_nominal_vs_calibrated", d.field), disc)

    # Auto-naming: every sample name contains the coordinate anchor
    for load in result.loads
        for (si, samp) in enumerate(load.samples)
            @test occursin("P$(lpad(si, 2,'0'))", samp.name)
        end
    end
end
```

- [ ] **Step 2: Run it**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/test/test_ingestion_core.jl
git commit -m "test(grouping): regression floor against SSRL 2026-04 Load 1 fixture"
```

---

## Task 10: Wire CLI `init` / `reingest` to call `scan_and_group!`

Replace the manifest-driven path in `cli_init_with_db!` and `reingest!` with calls to `scan_and_group!` when no manifest is present. The manifest path is retained as a legacy fallback (spec §9.5: "manifest support is retained only as a legacy fallback"). This means the CLI can ingest a beamtime directory without an `experiment.toml` manifest.

**Files:**
- Modify: `packages/HimalayaUI/src/cli.jl` (`cli_init_with_db!`, `_cli_init_inner!`, `reingest!`)
- Test: `packages/HimalayaUI/test/test_ingestion_core.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "cli_init_with_db! without manifest calls scan_and_group!" begin
    dir = mktempdir()
    data_dir     = joinpath(dir, "data")
    analysis_dir = joinpath(dir, "analysis")
    mkpath(data_dir); mkpath(analysis_dir)

    # Write an experiment.toml WITHOUT a manifest entry (or with manifest="")
    write(joinpath(dir, "experiment.toml"), """
    [experiment]
    name = "CLI Test"

    [layout]
    data_dir = "data"
    analysis_dir = "analysis"
    exposure_type = "simple"

    [files]
    integration = "{name}.dat"
    image = "{name}.tif"
    """)

    # Write 2 PRP/TIF pairs
    t0 = DateTime(2026, 4, 26, 23, 14, 8)
    for (i, h) in enumerate([58.9, 63.1])
        stem = "HA_$(i)_S$(lpad(i, 4,'0'))_0_001"
        write_prp(joinpath(data_dir, "$stem.prp");
            timestamp = Dates.format(t0 + Second((i-1)*19), "dd u yyyy HH:MM:SS"),
            beam_energy_ev = 9000.0, pipe_length_mm = 1700,
            detector = "Pilatus 1M", exposure_time = 15.0,
            horizontal_position_mm = h)
        write(joinpath(data_dir, "$stem.tif"), "fake tif")
    end
    write_setup_info(joinpath(analysis_dir, "setup_info_20260425_181705.txt"))

    db = fresh_db()
    exp_id = HimalayaUI.cli_init_with_db!(db, dir; analyze = false)

    # Exposures created via scan path
    exps = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, experiment_id FROM exposures WHERE experiment_id = ?", [exp_id]))
    @test length(exps) == 2

    # Samples created
    samps = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM samples WHERE experiment_id = ?", [exp_id]))
    @test length(samps) >= 1
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: FAIL (no exposures created — the current path requires a manifest).

- [ ] **Step 3: Modify `cli.jl`**

This is the one genuinely-conditional edit in the plan: the replacement target is a live function body, so line-verify before editing. As verified 2026-06-18, `_cli_init_inner!` (`cli.jl:41-118`) ends with this manifest-presence branch; the `else` body is a **single `println`** (`cli.jl:113-115`), preceded by the closing of the manifest `if`-block:

```julia
        println("Imported $sample_count samples and $exposure_count exposures from $(basename(manifest_path)).")
    else
        println("No manifest at $manifest_path — experiment registered without samples.")
    end

    println("Initialized experiment '$exp_name' (id=$exp_id) at $exp_dir")
    exp_id
end
```

Re-read the live block first (`grep -n 'No manifest at' packages/HimalayaUI/src/cli.jl`) and confirm the `else` body is still exactly that one `println`. Then replace **only** the `else` body, leaving the trailing `println(...)` / `exp_id` / `end` untouched:

```julia
        println("Imported $sample_count samples and $exposure_count exposures from $(basename(manifest_path)).")
    else
        # No manifest: run the automatic scan_and_group! ingest core (spec §9.5 —
        # manifest support retained only as a legacy fallback).
        println("No manifest at $manifest_path — running automatic directory scan...")
        scan_and_group!(db, exp_id, exp_dir;
            tif_pattern = cfg.image_pattern,
            dat_pattern = cfg.integration_pattern,
            analyze     = false)  # analysis runs outside the transaction in cli_init_with_db! (see below)
    end

    println("Initialized experiment '$exp_name' (id=$exp_id) at $exp_dir")
    exp_id
end
```

Variable names are verified live: `exp_id` (from `create_experiment!`, `cli.jl:63`), `exp_dir` (the function arg), `cfg` (`load_config` result, `cli.jl:53`), `manifest_path` (`cli.jl:58`).

> **`analyze = false` is correct:** `cli_init_with_db!` (`cli.jl:22-38`) calls `_analyze_experiment!(db, exp_id)` itself (`cli.jl:31`) after `_cli_init_inner!` returns and the transaction commits (only when `analyze=true` was passed to `cli_init_with_db!`). `scan_and_group!`'s own analysis step would double-analyze, so it is suppressed here.
>
> **Nested-transaction caveat (verify in Step 5's full-suite run):** `_cli_init_inner!` runs *inside* `SQLite.transaction(db)` (`cli.jl:24-25`), and `scan_and_group!` opens its own `SQLite.transaction(db)` blocks (Task 8). SQLite.jl's `transaction` does not nest as savepoints — a nested `transaction` call while one is already open can error or be a silent no-op depending on version. Two acceptable resolutions, pick one and note it in the commit: (a) call `scan_and_group!` here with the geometry/persist writes participating in the *outer* transaction by NOT opening inner `transaction` blocks when already in one (preferred — add an `in_transaction` guard kwarg to `scan_and_group!`), or (b) move the no-manifest `scan_and_group!` call out of `_cli_init_inner!` and into `cli_init_with_db!` *after* the outer transaction commits (alongside `_analyze_experiment!`). Confirm which the live SQLite.jl version requires before committing; the full-suite run in Step 5 will surface a nesting error if present.

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: PASS.

- [ ] **Step 5: Run the full backend suite once to catch regressions**

Run: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1; tail -40 /tmp/jl-test.out`
Expected: all green. Investigate any failure touching `create_experiment!`, `create_sample!`, or `create_exposure!` (likely callers that need the Phase-A signature sweep from Phase A Task 7/8).

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/cli.jl packages/HimalayaUI/test/test_ingestion_core.jl
git commit -m "feat(cli): cli_init_with_db! falls back to scan_and_group! when no manifest"
```

---

## Task 11: `ingest.jl` — `cheap_change_check` (Phase-C rescan-scheduler contract)

The Phase-C auto-rescan scheduler and the `POST /api/experiments/{id}/scan` route both call `HimalayaUI.cheap_change_check(db, experiment_id, root_dir)::Bool` and resolve it at call time via `isdefined(HimalayaUI, :cheap_change_check)` (Plan C `_rescan_tick!` `:1295`, `POST .../scan` `:1479`). **No plan defines it** — Plan C explicitly flags this as a Phase-B gap (`2026-06-18-ingestion-phase-c-scan-api-sse.md:1125-1131`, `:1763`): while it is absent, every scheduler tick and every manual scan is treated as "changed", forcing a full `scan_and_group!` on every tick and defeating the tiered-backoff "quiet directory" path (spec §9.4). This task closes that gap.

**Spec basis (§9.2, §9.4):** the rescan endpoint "runs a cheap change-check, then additive ingest … a no-change scan does nothing"; the scheduler tick "runs the cheap change-check first; only a changed directory triggers ingest". The spec does **not** prescribe a specific cheapness signal. **Signal chosen (minimal correct, no schema change):** compare the **count of matching image files in `root_dir`** against the **count of exposures already persisted** (`SELECT COUNT(*) FROM exposures WHERE experiment_id = ?`). The additive ingest dedups on `(experiment_id, filename)` (Phase A), so "more files on disk than exposures in the DB" is exactly "there are new files to ingest". This needs only a cheap `readdir`-style count (no `parse_prp`, no grouping, no per-file stat loop) and already-persisted state (the exposure count), so it requires **no new column**. A secondary newest-mtime signal is intentionally *not* used as a suppressor (mtime/`last_scanned_at` timezone parsing is fragile); the check is biased so that any uncertainty returns `true` — a false positive is a harmless extra scan (which is then a no-op via insert-only dedup), a false negative would silently drop new data and is unacceptable.

> **File choice — `ingest.jl`, not `grouping.jl`:** `cheap_change_check` reads the DB (`SQLite.DB` + a `COUNT(*)` query), which is orchestrator-level concern. `grouping.jl` is the pure, DB-free grouper (its docstring header: "no DB"); putting a DB query there would violate that layering. `ingest.jl` already owns the DB-touching orchestration (`scan_and_group!`) and is `include`d into the single `HimalayaUI` module, so `HimalayaUI.cheap_change_check` resolves module-wide exactly as Plan C's `isdefined` expects.

> **Optional future optimization (noted, not built here):** Phase A added `experiments.last_scanned_at` and the validated-2026-06-18 §9.4 note proposes `last_scan_tier` + `consecutive_empty_ticks` columns for the scheduler. A future cheaper signal could compare `max(mtime)` of the data dir against a stored `last_scanned_at` to detect *modified* (not just *added*) files without a count. That is deliberately out of scope here — the count signal is sufficient for the additive-ingest contract (modified-but-not-added files do not change grouping), and adding mtime comparison now would risk false negatives from timezone skew.

**Files:**
- Modify: `packages/HimalayaUI/src/ingest.jl` (append `cheap_change_check`)
- Test: `packages/HimalayaUI/test/test_ingestion_core.jl`

- [ ] **Step 1: Write the failing test**

Append to the `@testset "ingestion core (Phase B)"` block. Reuses the `write_prp` / `write_setup_info` / `fresh_db` helpers from Task 0 and the same fixture shape as the Task 8 test.

```julia
@testset "cheap_change_check" begin
    dir = mktempdir()
    data_dir     = joinpath(dir, "data")
    analysis_dir = joinpath(dir, "analysis")
    mkpath(data_dir); mkpath(analysis_dir)

    t0 = DateTime(2026, 4, 26, 23, 14, 8)
    # Write 2 TIF/PRP pairs
    for (i, h) in enumerate([58.9, 63.1])
        stem = "HA_$(i)_S$(lpad(i, 4,'0'))_0_001"
        write_prp(joinpath(data_dir, "$stem.prp");
            timestamp = Dates.format(t0 + Second((i-1)*19), "dd u yyyy HH:MM:SS"),
            beam_energy_ev = 9000.0, pipe_length_mm = 1700,
            detector = "Pilatus 1M", exposure_time = 15.0,
            horizontal_position_mm = h)
        write(joinpath(data_dir, "$stem.tif"), "fake tif")
    end
    write_setup_info(joinpath(analysis_dir, "setup_info_20260425_181705.txt"))

    db = fresh_db()
    exp_id = HimalayaUI.create_experiment!(db;
        name = "cheap-check", path = dir,
        data_dir = data_dir, analysis_dir = analysis_dir)

    # Before any ingest: 2 files on disk, 0 exposures in DB → changed
    @test HimalayaUI.cheap_change_check(db, exp_id, dir) == true

    # After ingest: 2 files, 2 exposures → unchanged (true no-op tick)
    HimalayaUI.scan_and_group!(db, exp_id, dir; analyze = false)
    @test HimalayaUI.cheap_change_check(db, exp_id, dir) == false

    # Add a new TIF/PRP pair on disk (not yet ingested) → changed again
    new_stem = "HA_3_S0003_0_001"
    write_prp(joinpath(data_dir, "$new_stem.prp");
        timestamp = Dates.format(t0 + Second(600), "dd u yyyy HH:MM:SS"),
        beam_energy_ev = 9000.0, pipe_length_mm = 1700,
        detector = "Pilatus 1M", exposure_time = 15.0,
        horizontal_position_mm = 67.3)
    write(joinpath(data_dir, "$new_stem.tif"), "fake tif")
    @test HimalayaUI.cheap_change_check(db, exp_id, dir) == true

    # Re-scan picks up the new file; check goes quiet again
    HimalayaUI.scan_and_group!(db, exp_id, dir; analyze = false)
    @test HimalayaUI.cheap_change_check(db, exp_id, dir) == false

    # Defensive: a non-existent data dir is treated as "no change" (nothing to ingest),
    # never an error (the scheduler tick must not crash on a vanished volume).
    bad_db = fresh_db()
    bad_id = HimalayaUI.create_experiment!(bad_db;
        name = "missing-dir", path = dir,
        data_dir = joinpath(dir, "does_not_exist"),
        analysis_dir = analysis_dir)
    @test HimalayaUI.cheap_change_check(bad_db, bad_id, dir) == false
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: FAIL (`UndefVarError: cheap_change_check not defined`).

- [ ] **Step 3: Append `cheap_change_check` to `ingest.jl`**

```julia
"""
    cheap_change_check(db, experiment_id, root_dir) -> Bool

Cheap "has the directory changed since the last ingest?" probe for the Phase-C
auto-rescan scheduler and the `POST /api/experiments/{id}/scan` route (both
resolve this by name via `isdefined(HimalayaUI, :cheap_change_check)`).

Returns `true` when the data directory appears to hold files not yet persisted
(so a `scan_and_group!` is warranted), `false` when it looks unchanged (the
scheduler tick can stay quiet / back off — spec §9.4).

**Cheapness contract:** this does NOT parse PRP files or run grouping. It counts
matching image files in the experiment's `data_dir` (a single `readdir` + suffix
filter) and compares against `COUNT(*)` of already-persisted exposures. Additive
ingest dedups on `(experiment_id, filename)`, so "more files on disk than rows in
the DB" is exactly "there is new data to ingest".

**Bias:** any ambiguity returns `true` (an extra scan is a harmless no-op via the
insert-only dedup; a missed scan would silently drop data). The only `false`
shortcuts are (a) a vanished/unreadable data dir — nothing to ingest, and a
scheduler tick must never crash on a missing volume — and (b) on-disk count ≤
persisted count.

`root_dir` is accepted to match the Phase-C call contract
(`cheap_change_check(db, experiment_id, root_dir)`); the authoritative scan root
is the experiment's stored `data_dir`, so `root_dir` is currently advisory only.
"""
function cheap_change_check(
    db::SQLite.DB,
    experiment_id::Int,
    root_dir::AbstractString;
    image_pattern::String = "{name}.tif",
)::Bool
    # Resolve the experiment's data_dir (the authoritative scan root).
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT data_dir FROM experiments WHERE id = ?", [experiment_id]))
    isempty(rows) && return true   # unknown experiment: let the caller's scan surface the error
    data_dir = String(first(rows).data_dir)

    isdir(data_dir) || return false  # vanished/unreadable dir: nothing to ingest, never crash

    # Cheap on-disk count: number of files whose name matches the image suffix.
    # image_pattern is "{name}<suffix>"; everything after "{name}" is the literal suffix.
    suffix = replace(image_pattern, "{name}" => "")
    on_disk = try
        count(f -> endswith(f, suffix), readdir(data_dir))
    catch
        return true   # readdir failed for any reason: assume changed (safe direction)
    end

    # Persisted count.
    persisted = Int(first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM exposures WHERE experiment_id = ?", [experiment_id]))).n)

    return on_disk > persisted
end
```

> **Pattern note:** the suffix-strip mirrors how Phase B's `scan_directory` enumerates (it passes the same `image_pattern` to `resolve_files`). The default `"{name}.tif"` matches `scan_and_group!`'s default; the scheduler may pass the experiment's configured `image_pattern` if it diverges (Plan C currently calls the 3-arg positional form, so the default applies — adequate for the SSRL `.tif` convention). Keeping `image_pattern` as a defaulted kwarg means the Phase-C 3-arg call site compiles unchanged.

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/ingest.jl packages/HimalayaUI/test/test_ingestion_core.jl
git commit -m "feat(ingest): cheap_change_check for the Phase-C rescan scheduler"
```

---

## Task 12: `grouping.jl` — `derive_sample_flags` (read-time merge/split flag derivation)

Expose the **pure** `derive_sample_flags(load_rows) → Dict{Int, GroupingFlag}` keyed by `sample_id` (spec §8.8, §9.1). It is a **read-time derivation over the grouped rows** — NO DB writes, NOT a stored column, NOT part of `scan_and_group!`'s persistence. Phase D's `get_loads_rollup` (`2026-06-18-ingestion-phase-d-structural-events.md:456`) is the cross-plan consumer: it builds the nested rows from the DB, then calls `derive_sample_flags` over them and suppresses any flag with a non-undone `grouping_flag_dismissed` event (spec §9.2, §9.3).

**Input contract — the `get_loads_rollup` shape, verified 2026-06-18** (`…phase-d-structural-events.md:443-482`). `derive_sample_flags` consumes exactly what that function returns: a `Vector{<:NamedTuple}` of loads, each `(load_id, load_index, start_time, end_time, frame_count, samples)`; each sample `(sample_id, name, slot_index, grouping_source, name_source, merged_into_id, exposures)`; each exposure `(id, filename, horizontal_position, timestamp)` (the §8.8 `LoadExposure` leaf — keyed `id`, NOT `exposure_id`; no `frame_no`/`status`). `derive_sample_flags` reads only `horizontal_position` + `filename` per exposure (and the sample's `sample_id`/`name`), so it is agnostic to any extra fields the roll-up may also carry. So:
- the **per-exposure horizontal-position field is `horizontal_position`** (the persisted/rollup name — *not* the pre-persistence `GroupedExposure.horizontal_position_mm` from Task 7), nullable;
- the **per-sample label is `name`** (composed `<label> (SNNPMM)`, with a per-load-varying coordinate, so it is *not* a cross-load match key); the **cross-load shared label** is recovered from each sample's exposures' `filename` via the existing `_label_from_stem` (Task 7);
- `merge_with_sample_id` / `merge_with_label` and `split_at_index` / `jump_from` / `jump_to` come from those fields.

**Julia representation of `GroupingFlag` (serializes to the §8.8 contract JSON):** a tiny struct union, returned only for flagged samples (an unflagged sample is *absent* from the Dict — the contract's `null`):

```julia
struct MergeFlag
    merge_with_sample_id ::Int
    merge_with_label     ::String
end
struct SplitFlag
    split_at_index ::Int
    jump_from      ::Float64
    jump_to        ::Float64
end
const GroupingFlag = Union{MergeFlag, SplitFlag}
```

**How Phase D serializes it** (cross-plan note, not built here): in `get_loads_rollup`'s per-sample JSON, the `flag` key is `nothing` when the sample id is absent from the Dict (→ JSON `null`), `Dict("kind"=>"merge", "merge_with_sample_id"=>f.merge_with_sample_id, "merge_with_label"=>f.merge_with_label)` for a `MergeFlag`, and `Dict("kind"=>"split", "split_at_index"=>f.split_at_index, "jump_from"=>f.jump_from, "jump_to"=>f.jump_to)` for a `SplitFlag` — matching the §8.8 TS `GroupingFlag` discriminated union verbatim. (A sample may have at most one flag in this Dict; if both a merge candidate and a split are detected for the same sample, **split wins** — an internal split must be resolved before a cross-load merge is meaningful.)

**Files:**
- Modify: `packages/HimalayaUI/src/grouping.jl` (append `MergeFlag`/`SplitFlag`/`GroupingFlag` + `derive_sample_flags`)
- Test: `packages/HimalayaUI/test/test_ingestion_core.jl`

- [ ] **Step 1: Write the failing test**

The fixture builds `get_loads_rollup`-shaped NamedTuples directly (no DB) and exercises **both** a cross-load merge candidate (the same filename label `HA_85` recurs in Load 1 and Load 2) **and** an intra-sample split jump (a sample whose two exposures sit ~4 mm apart — beyond the local within-burst jitter — including the pure single-frame round-robin case the prior `_cluster_slots` deferred to "flag for review"). Append to the `@testset "ingestion core (Phase B)"` block:

```julia
@testset "derive_sample_flags" begin
    # Helper to build a rollup-shaped exposure NamedTuple (mirrors get_loads_rollup).
    exp(id, fname, hpos) = (id = id, filename = fname,
        horizontal_position = hpos, timestamp = nothing)
    smp(sid, name, slot, exps) = (sample_id = sid, name = name, slot_index = slot,
        grouping_source = "auto_position", name_source = "auto",
        merged_into_id = nothing, exposures = exps)
    load(lid, lidx, smps) = (load_id = lid, load_index = lidx,
        start_time = nothing, end_time = nothing,
        frame_count = sum(length(s.exposures) for s in smps), samples = smps)

    # --- Cross-load MERGE candidate: label "HA_85" recurs across two loads ---
    # Load 1 slot 1: HA_85 burst (clean, ~0.1 mm jitter) → merge candidate (recurs in L2)
    # Load 2 slot 1: HA_85 burst (clean)                  → merge candidate (recurs in L1)
    # Load 1 slot 2: HA_90 burst (no recurrence)          → no flag
    s1 = smp(101, "HA85 (S01P01)", 1,
             [exp(1, "HA_85_S2001_0_001", 74.80), exp(2, "HA_85_S2002_0_002", 74.83)])
    s2 = smp(102, "HA90 (S01P02)", 2,
             [exp(3, "HA_90_S2003_0_001", 70.85), exp(4, "HA_90_S2004_0_002", 70.88)])
    s3 = smp(201, "HA85 (S02P01)", 1,
             [exp(5, "HA_85_S2101_0_001", 74.75), exp(6, "HA_85_S2102_0_002", 74.79)])

    rollup_merge = [load(1, 1, [s1, s2]), load(2, 1, [s3])]
    flags = HimalayaUI.derive_sample_flags(rollup_merge)

    # s1 and s3 are merge candidates of each other (shared label "HA_85")
    @test haskey(flags, 101)
    @test flags[101] isa HimalayaUI.MergeFlag
    @test flags[101].merge_with_sample_id == 201
    @test flags[101].merge_with_label == "HA_85"
    @test haskey(flags, 201)
    @test flags[201].merge_with_sample_id == 101
    # s2 (label HA_90, no recurrence) is NOT flagged → absent from the Dict (contract null)
    @test !haskey(flags, 102)

    # --- Intra-sample SPLIT: one sample spans a ~4 mm position jump ---
    # Load 1 has only this sample, so there is no cross-load recurrence;
    # within the sample the position jumps 74.80 → 67.20 (≫ local jitter).
    split_sample = smp(301, "HA77 (S01P01)", 1,
        [exp(10, "HA_77_S3001_0_001", 74.80),
         exp(11, "HA_77_S3002_0_002", 74.85),   # still ~slot center
         exp(12, "HA_77_S3003_0_003", 67.20),   # JUMP here (index 3)
         exp(13, "HA_77_S3004_0_004", 67.25)])
    rollup_split = [load(3, 1, [split_sample])]
    sflags = HimalayaUI.derive_sample_flags(rollup_split)
    @test haskey(sflags, 301)
    @test sflags[301] isa HimalayaUI.SplitFlag
    @test sflags[301].split_at_index == 3            # the exposure index where the jump occurs
    @test sflags[301].jump_from ≈ 74.85
    @test sflags[301].jump_to   ≈ 67.20

    # --- Pure single-frame round-robin (the prior Plan-B deferral, now a concrete split) ---
    # Each "exposure" is a different slot position visited in one round; _cluster_slots
    # under-split it to ONE sample (KNOWN LIMITATION). derive_sample_flags flags the
    # first position jump as a split suggestion (spec §5: "single-frame … else flag").
    rr_sample = smp(401, "JC (S01P01)", 1,
        [exp(20, "JC_C01_S4001_0_001", 39.5),
         exp(21, "JC_C02_S4002_0_001", 71.0),    # jump → split at index 2
         exp(22, "JC_C03_S4003_0_001", 87.0)])
    rr_flags = HimalayaUI.derive_sample_flags([load(4, 1, [rr_sample])])
    @test haskey(rr_flags, 401)
    @test rr_flags[401] isa HimalayaUI.SplitFlag
    @test rr_flags[401].split_at_index == 2

    # --- No flags for a clean directory (one load, one clean burst) ---
    clean = smp(501, "HA60 (S01P01)", 1,
        [exp(30, "HA_60_S5001_0_001", 58.9), exp(31, "HA_60_S5002_0_002", 58.95)])
    @test isempty(HimalayaUI.derive_sample_flags([load(5, 1, [clean])]))
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: FAIL (`UndefVarError: derive_sample_flags not defined`).

- [ ] **Step 3: Append `derive_sample_flags` to `grouping.jl`**

```julia
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
    derive_sample_flags(load_rows; slot_k = 5.0) -> Dict{Int, GroupingFlag}

PURE read-time derivation of per-sample merge/split suggestions over the
`get_loads_rollup` rows (spec §8.8). No DB access; returns a Dict keyed by
`sample_id`, with a sample **absent** from the Dict meaning "no flag" (the
contract's JSON `null`).

`load_rows` is the `get_loads_rollup` shape (verified 2026-06-18):
`Vector` of loads `(load_id, load_index, …, samples)`; each sample
`(sample_id, name, slot_index, …, exposures)`; each exposure
`(id, filename, horizontal_position, timestamp)` (the §8.8 leaf).

Two suggestion kinds (a sample gets at most one; **split wins** if both apply):

1. **Split** — within one sample, the per-exposure `horizontal_position` jumps
   beyond the *local within-burst jitter* (gap-relative per spec §5, derived per
   load, NOT a fixed 4 mm). The jitter tolerance mirrors `_cluster_slots`'s rule:
   `slot_k × median(|consecutive deltas|)`, with the same median-near-zero
   fallback (learn from non-zero deltas; `Inf` when there are none). This is also
   where the **pure single-frame round-robin** case the grouper deferred (Task 6
   KNOWN LIMITATION — `_cluster_slots` under-splits it to one slot) becomes a
   concrete split suggestion: each exposure sits at a different slot position, so
   the first position gap exceeding the (per-load) tolerance is flagged.

2. **Merge** — a sample's filename label (via `_label_from_stem` over its
   exposures' `filename`s) recurs as the label of a sample in *another* load. The
   grouper never auto-merges cross-load (spec §5), so this is surfaced as a
   suggestion pointing at the other sample (`merge_with_sample_id` +
   `merge_with_label`). When a label recurs in more than two loads, each flagged
   sample points at the *first other* sample sharing that label (lowest
   `sample_id`); the UI walks the chain.

`slot_k = 5.0` matches `_cluster_slots`.
"""
function derive_sample_flags(load_rows; slot_k::Float64 = 5.0)::Dict{Int, GroupingFlag}
    flags = Dict{Int, GroupingFlag}()

    # ---- Per-load jitter tolerance, mirroring _cluster_slots (spec §5) -----
    # Returns the slot-spacing tolerance for one load from its consecutive-frame
    # position deltas. Inf when there is no movement at all (one slot).
    function _load_jitter_tol(positions::Vector{Union{Float64, Missing}})
        deltas = Float64[]
        for i in 2:length(positions)
            if !ismissing(positions[i]) && !ismissing(positions[i-1])
                push!(deltas, abs(positions[i] - positions[i-1]))
            end
        end
        isempty(deltas) && return Inf
        med = Statistics.median(deltas)
        if med < 1e-6
            # Fallback regime: nonzero deltas ARE slot spacing → tolerance below them
            # (÷ slot_k), not above (see Task 6 fallback note; corrected 2026-06-19).
            nonzero = filter(d -> d > 1e-6, deltas)
            return isempty(nonzero) ? Inf : Statistics.median(nonzero) / slot_k
        end
        return med * slot_k
    end

    # =====================================================================
    # 1. SPLIT suggestions — per sample, per load (uses the load-local jitter
    #    computed over *all* of that load's exposures in row order).
    # =====================================================================
    for ld in load_rows
        # Flatten the load's exposures in (slot, then row) order to learn the
        # load-local position-delta distribution — the same population
        # _cluster_slots used during ingest.
        load_positions = Union{Float64, Missing}[]
        for sm in ld.samples, ex in sm.exposures
            push!(load_positions,
                ex.horizontal_position === nothing ? missing : Float64(ex.horizontal_position))
        end
        tol = _load_jitter_tol(load_positions)

        for sm in ld.samples
            xs = sm.exposures
            length(xs) < 2 && continue
            for i in 2:length(xs)
                p_prev = xs[i-1].horizontal_position
                p_curr = xs[i].horizontal_position
                (p_prev === nothing || p_curr === nothing) && continue
                if abs(Float64(p_curr) - Float64(p_prev)) > tol
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
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_ingestion_core.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/grouping.jl packages/HimalayaUI/test/test_ingestion_core.jl
git commit -m "feat(grouping): derive_sample_flags pure read-time merge/split derivation"
```

> **Cross-plan consumer (Phase D):** `get_loads_rollup` (`2026-06-18-ingestion-phase-d-structural-events.md:456`) calls `derive_sample_flags(rollup)` over the rows it just read, attaches each sample's flag (serialized per the §8.8 union, `null` when absent), and **suppresses** any flag matching a non-undone `grouping_flag_dismissed` event (spec §9.2/§9.3) before returning JSON. That suppression is Phase D's responsibility — `derive_sample_flags` here is the pure producer and has no knowledge of dismiss events.

---

## Self-Review

### Spec coverage (§5, §6, §9.1)

| Spec requirement | Task |
|---|---|
| §9.1: `prp.jl` / `parse_prp(path) → NamedTuple`, defensive | Task 1 ✓ |
| §9.1: `geometry.jl` / `derive_geometry(prp_paths, setup_files)` | Tasks 2–3 ✓ |
| §6: Detector → pitch lookup, unknown → missing + flag | Task 2 ✓ |
| §6: Setup file `Beam center is at (X, Y)` + `Mean distance` parse | Task 2 ✓ |
| §6: Geometry authority — setup file `Mean distance` > PRP `Pipe length` | Task 3 ✓ |
| §6: Per-field `*_source ∈ {prp, setup, human, default}` | Task 3 ✓ |
| §6: Multi-setup discrepancy detection (beam energy, pipe length variance) | Task 3 ✓ |
| §6: Latest setup file by filename sort | Task 3 ✓ |
| §9.1: `grouping.jl` / `scan_directory` reusing `resolve_files` | Task 4 ✓ |
| §9.1: dispatch gap resolved (minimal ExperimentConfig for dispatch) | Task 4 ✓ |
| §5: Time-gap load segmentation (gap-relative, not absolute) | Task 5 ✓ |
| §5: Unimodal fallback (no clear separation → one load + flag) | Task 5 ✓ (flag now propagated into `GroupedLoad.flag` + `GroupingResult.discrepancies` in Task 7) |
| §5: Slot clustering keyed on "gap ≫ median consecutive-frame delta" | Task 6 ✓ |
| §5: Median-near-zero fallback (multi-frame bursts → learn gap from non-zero deltas) | Task 6 ✓ |
| §5: Pure single-frame round-robin (no burst anywhere) → flag for review | Task 6 under-splits to one slot during ingest (KNOWN LIMITATION), but Task 12's `derive_sample_flags` converts that into a concrete `SplitFlag` at read time — so the case **is** flagged for the human, just at the rollup layer rather than during clustering |
| §4/§10: persist `scan_id` / `frame_no` (parsed from filename stem) | Tasks 7–8 ✓ (`_parse_scan_frame`) |
| §4: persist `image_path` (image route + `prewarm_thumbnails!` require it) | Task 8 ✓ (`image_path = ge.tif_path`) |
| §5: Stepping axis = motor with dominant inter-burst variance (horizontal_position_mm) | Task 6 ✓ (detected from data, not hard-coded) |
| §5: `group_into_samples` top-level grouper | Task 7 ✓ |
| §5: Naming rule `<label> (SNNPMM)`, `_label_from_stem` extraction | Task 7 ✓ |
| §9.1: `scan_and_group!` transactional orchestrator, `_DB_WRITE_LOCK` | Task 8 ✓ |
| §9.1: insert-only discipline (dedup by `(experiment_id, filename)`) | Task 8 ✓ |
| §4: never-clobber for `name_source='human'` geometry fields | Task 8 (`_update_geometry_if_not_human!`) ✓ |
| §9.1: `analyze_exposure!` called outside main transaction | Task 8 ✓ |
| §9.5: CLI `init` rewrite to call `scan_and_group!` core | Task 10 ✓ |
| §9.2/§9.4: cheap change-check for rescan endpoint + scheduler (`cheap_change_check`) | Task 11 ✓ (closes the Phase-C cross-plan gap) |
| §8.8/§9.1: `grouping.jl` exposes pure `derive_sample_flags(load_rows) → Dict{sample_id, GroupingFlag}` (merge + split, no DB) | Task 12 ✓ (consumed cross-plan by Phase D's `get_loads_rollup`, `…phase-d-structural-events.md:456`, which suppresses dismissed flags) |
| §5: Pure single-frame round-robin → concrete split suggestion (the Task-6 deferral) | Task 12 ✓ (now a `SplitFlag`; the read-time derivation flags what `_cluster_slots` under-split) |
| §11: regression-floor assertions, real fixture parameters | Task 9 ✓ |

### Deferred to later plans (not gaps)
- Phase C: `POST /api/experiments`, `POST /api/experiments/{id}/scan`, `GET /api/experiments/{id}/loads`, `broadcast_progress!`.
- Phase D: `exposure_moved`, `sample_renamed`, `sample_created`, `sample_split`, `grouping_flag_dismissed` event kinds (no `sample_merged`).
- Phase E: frontend `/experiments/:id` shell, `ExperimentsHomePage`, `GroupingReviewPage`, all new React components.
- Counter-reset-to-1 as a *corroborating* hint for cross-load identity (spec §5): demoted by spec; implementation deferred to Phase D when the grouping-review UI is available.
- **Pure single-frame round-robin clustering (spec §5 "single-frame acquisitions … else flag"):** a load with *no* multi-frame burst anywhere (each slot exactly one frame, slots revisited round-robin) cannot have its slot tolerance learned from within-burst jitter — the jitter and the slot-spacing are the same delta population, so `_cluster_slots` under-splits to a single slot. The spec's prescribed behavior is to **flag for human review**, which requires the Phase-D grouping-review UI to surface and resolve. For Phase B this case falls through to one slot (documented in the `_cluster_slots` docstring KNOWN LIMITATION). The Phase-B test for the fallback therefore exercises the case the algorithm *does* handle (multi-frame bursts → median-near-zero → learn-from-non-zero-deltas), not the unhandled round-robin case.
- CLI `reingest!` rewrite to `scan_and_group!` (the current `reingest!` requires `_reingest_inner!` which needs a manifest; the path to `scan_and_group!` is Task 10's `cli_init_with_db!` case; the full `reingest!` swap is done in Phase C alongside the `POST /api/experiments/{id}/scan` route).

### Placeholder scan
Every code step contains complete, runnable Julia code. No "TBD", no "similar to Task N", no "add error handling" without the code. The four `@warn` calls in `parse_prp` / `parse_setup_info` / `scan_and_group!` are the actual error-handling. Task 10's "read the live block first" note is an accuracy guard (it targets a live function body), not a placeholder — both the current live `else` body and the full replacement are given verbatim. The Task 10 nested-transaction caveat names two concrete resolutions (it is a genuine SQLite.jl-version-dependent decision, surfaced for the executor, not an omission).

### Type/name consistency
- `parse_prp(path) → NamedTuple` — used identically in its definition (Task 1), `scan_directory` (Task 4), and `group_into_samples` (Task 7 via `ExposureMeta.prp`).
- `derive_geometry(prp_paths, setup_files) → (geometry::NamedTuple, discrepancies::Vector{GeometryDiscrepancy})` — consistent across Tasks 3, 8, 9.
- `scan_directory(data_dir, analysis_dir) → Vector{ExposureMeta}` — consistent across Tasks 4, 8, 9.
- `group_into_samples(metas) → GroupingResult` — consistent across Tasks 7, 8, 9.
- `scan_and_group!(db, exp_id, dir; ...) → NamedTuple` — consistent across Task 8, Task 10.
- `cheap_change_check(db::SQLite.DB, experiment_id::Int, root_dir::AbstractString)::Bool` (Task 11) — **matches the Phase-C call contract exactly**: Plan C calls `HimalayaUI.cheap_change_check(db, experiment_id, root_dir)` (3-arg positional) resolved via `isdefined(HimalayaUI, :cheap_change_check)` at `_rescan_tick!` (`2026-06-18-ingestion-phase-c-scan-api-sse.md:1295-1296`) and the rescan route (`:1479`). The extra `image_pattern` kwarg is defaulted, so the 3-arg positional call site compiles unchanged.
- `_segment_loads_with_flag(metas) → (loads, flag)` / `_segment_loads(metas)` — both defined Task 5; `_cluster_slots(load_metas)` defined Task 6; all called in Task 7.
- `GroupedExposure`, `GroupedSample`, `GroupedLoad`, `GroupingResult` defined Task 7, used in Task 8.
- `MergeFlag` / `SplitFlag` / `const GroupingFlag = Union{MergeFlag, SplitFlag}` and `derive_sample_flags(load_rows; slot_k) → Dict{Int, GroupingFlag}` defined Task 12. Input is the **`get_loads_rollup` row shape** (Phase D, `…phase-d-structural-events.md:443-482`) — per-exposure position field is `horizontal_position` (rollup name), per-sample label recovered via `_label_from_stem` (Task 7) over exposure `filename`s. Phase D's `get_loads_rollup` is the sole consumer and owns the §8.8 JSON serialization + dismissed-flag suppression.
- `GeometryDiscrepancy` defined Task 3, referenced in Task 9 test.
- `ExposureMeta` defined Task 4, used throughout.

### Known accuracy risks
- `ExperimentConfig` positional constructor in `_minimal_scan_config` (Task 4): `ExperimentConfig` is a plain `struct` with **no keyword constructor** (`config.jl:21`, verified 2026-06-18), so only the 23-field positional form is valid — the Task 4 call matches the live field order exactly. Re-verify against `config.jl:21-50` before editing if the struct has since changed.
- `Dates.value(ts_curr - ts_prev)` returns milliseconds for `DateTime` differences in Julia; the division by 1000.0 in Task 5 is the correct ms→s conversion (confirmed by `Dates.value(Second(1)) == 1000`).
- `scan_and_group!` calls `create_experiment!` with `path=dir` — confirm the Phase-A `create_experiment!` signature takes `path::String` (it does, at `db.jl:1810`).
- Task 10's `scan_and_group!` call inside `_cli_init_inner!` passes `analyze=false` to avoid double-analysis; confirm the outer `cli_init_with_db!` still calls `_analyze_experiment!` after the `_cli_init_inner!` transaction completes (it does, at `cli.jl:31`). See also the nested-transaction caveat in Task 10 Step 3.

### Review applied (2026-06-18)

A review pass against live source folded in the following fixes (all re-verified against `db.jl`/`config.jl`/`cli.jl`/`image.jl` and the corrected Phase-A plan before writing):

- **[P0] `[][]` empty-return literal** (BoundsError): `Vector{ExposureMeta}[][]` is an empty `Vector{Vector{…}}` immediately `getindex`-ed → runtime BoundsError. Both occurrences fixed — `_segment_loads_with_flag` empty branch now returns `(Vector{ExposureMeta}[], :ok)`, `_cluster_slots` empty branch returns `Vector{ExposureMeta}[]`.
- **[P0] `scan_and_group!` rescan idempotency**: the per-load `INSERT INTO loads` was unconditional; on rescan it minted a fresh `load_id`, re-keyed the sample dedup, and duplicated samples + orphaned the load row (the `loads` table has only a non-unique index — Phase A Task 2). Fixed with a `SELECT id FROM loads WHERE experiment_id=? AND load_index=?` reuse-or-insert (mirrors the sample dedup). Task 8 test now also asserts loads and samples row counts are unchanged after the second scan.
- **[P1] `image_path` never persisted**: added `image_path = ge.tif_path` to the `create_exposure!` call (the image route + `prewarm_thumbnails!` both filter `WHERE image_path IS NOT NULL`, and the legacy manifest path always sets it). Task 8 asserts non-NULL `.tif` `image_path`.
- **[P1] `scan_id` / `frame_no` never persisted**: added `_parse_scan_frame(stem)` (parses the `_S<digits>_` token + trailing frame index from the filename stem), two new `GroupedExposure` fields, and the `scan_id=`/`frame_no=` kwargs on the `create_exposure!` call. Task 8 asserts both are persisted.
- **[P1] Segmentation review-flag discarded**: `group_into_samples` now binds `seg_flag` from `_segment_loads_with_flag` and propagates it into `GroupedLoad.flag` and (for `:unimodal_fallback`) into `GroupingResult.discrepancies`, instead of hard-coding `:ok`.
- **[P1] Single-frame fallback test contradiction**: the round-robin fixture could not be split by the documented algorithm (jitter == spacing population), so the `>= 2` assertion was unreachable. Rewrote the fallback test to exercise the case the `med_delta < 1e-6` branch actually handles (multi-frame bursts → learn tolerance from non-zero deltas, asserting exactly 3 slots), aligned the docstring + inline comment, and documented the genuinely-unhandled pure-single-frame round-robin case as a Phase-D deferral (flag-for-review per spec).
- **[P2] `_minimal_scan_config` q_units** `"A^-1"` → `"A-1"` (codebase default, `config.jl:112`); replaced the misleading `@kwdef`-style keyword-constructor fallback note — `ExperimentConfig` is a plain `struct` (no keyword constructor; would `MethodError`) with a verified 23-field positional order.
- **[P2] `write_prp` dead kwargs**: dropped the unused `scan_id`/`frame_no` kwargs from `write_prp` and its two callers (scan/frame now come from the filename stem, consistent with `_parse_scan_frame`).
- **[P2] `runtests.jl` anchor**: corrected the stale "test_validate.jl is the last entry" claim (the live last entries are `test_migrate_toml.jl` then `test_spa_fallback.jl`).
- **[P2] Task 10 manifest else-branch**: replaced the elided/hypothetical edit with the exact live `else` body + the exact replacement, verified variable names, and surfaced the nested-`SQLite.transaction` caveat (the one genuinely-conditional edit).
