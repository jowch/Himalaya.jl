# Beam-center Calibration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let a scientist declare the beam center + detector pixel size in `experiment.toml` so the Focus page renders physically-correct, beam-centered q-rings (instead of today's centered/decorative fallback).

**Architecture:** Three new optional `[beamline]` fields (`beam_center_x`, `beam_center_y`, `pixel_size_um`) round-trip through `ExperimentConfig` and are surfaced (with the existing `energy_kev`/`flight_path_m`) on the experiment JSON — parsed from the config blob, no DB migration (the `q_units` precedent). The image route emits the raw detector pixel dimensions as response headers (the displayed PNG is downscaled, so the bitmap size is not raw). A pure frontend helper assembles a `DetectorCalibration` and feeds the **already-complete** geometry engine (`detectorGeometry.ts` + `DetectorRings`). Incomplete calibration → unchanged centered fallback.

**Tech Stack:** Julia (HimalayaUI backend, SQLite, Oxygen.jl, TOML, stdlib Test), React/TypeScript (Vite, Vitest, TanStack Query).

**Spec:** [docs/superpowers/specs/2026-06-14-beam-center-calibration-design.md](../specs/2026-06-14-beam-center-calibration-design.md) — read it first; it carries the design rationale and the three verification findings (canvas-not-`<img>`, full-route resize, landscape rotation).

---

## File Structure

**Backend (Julia):**
- Modify: `packages/HimalayaUI/src/config.jl` — 3 struct fields + `_build_config` reads + `config_to_toml` writes
- Modify: `packages/HimalayaUI/src/routes_experiments.jl` — add `_beamline_from_config`; `_q_units_from_config` delegates to it; `_experiment_row_to_json` merges 3 new fields
- Modify: `packages/HimalayaUI/src/routes_exposures.jl` — image route emits `X-Image-Width`/`X-Image-Height` (full branch only)
- Test: `packages/HimalayaUI/test/test_config.jl`, `test_routes_experiments.jl`, `test_routes_image.jl`

**Frontend (TypeScript):**
- Modify: `packages/HimalayaUI/frontend/src/api.ts` — `Experiment` type: 5 fields
- Modify: `packages/HimalayaUI/frontend/src/print/pages/focusAdapters.ts` — new `buildDetectorCalibration` pure helper
- Modify: `packages/HimalayaUI/frontend/src/print/detector/DetectorImage.tsx` — `onRawSize`/`onOrient` callbacks
- Modify: `packages/HimalayaUI/frontend/src/print/components/DetectorPanel.tsx` — forward callbacks + pass `orient` to `DetectorRings`
- Modify: `packages/HimalayaUI/frontend/src/print/pages/FocusPage.tsx` — state + calibration + `imageAspect` wiring
- Test: `frontend/test/print-pages/focusAdapters.test.ts`, `frontend/test/print-detector/DetectorImage.test.tsx`, `frontend/test/print-components/DetectorPanel.test.tsx`

**Docs/template:**
- Modify: `packages/HimalayaUI/configs/simple.toml`, `docs/experiment-config.md`

**Dependency order:** Task 1 → 2 → 3 (backend, independent of frontend) ‖ Task 4 → 5 → 6 → 7 → 8 (frontend chain) → 9 (docs) → 10 (full gates + live verify). Backend and frontend chains are independent until Task 10.

### Running tests (read before starting)

- **Julia, single file (fast iteration):**
  - `test_config.jl` is self-contained: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_config.jl`
  - Route/image tests need the `with_test_server` helper from `test_http.jl`. Run from the package dir:
    ```bash
    cd packages/HimalayaUI && julia --project=. -e 'using Test, HimalayaUI, SQLite, HTTP, JSON3, DBInterface, Tables, FileIO, ImageCore, TiffImages; include("test/test_http.jl"); include("test/test_routes_experiments.jl")'
    ```
    (swap in `test/test_routes_image.jl` for Task 3). `test_http.jl` runs its own testset too — that's fine.
- **Julia, full suite (authoritative gate, ~5–10 min — Task 10 only, capture once):**
  ```bash
  julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
  grep -E "Test Summary|did not pass|fail" /tmp/jl-test.out ; tail -50 /tmp/jl-test.out
  ```
  Do NOT re-run the full suite with different greps — it rebuilds fixtures each time.
- **Frontend, single file:** from `packages/HimalayaUI/frontend/`: `npm test -- <substring>` (one-shot; e.g. `npm test -- focusAdapters`).
- **Frontend gates (Task 10):** `npm run build` (tsc --noEmit + vite), `npm run lint:design`.

---

## Task 1: Config — `[beamline]` beam center + pixel size round-trip

**Files:**
- Modify: `packages/HimalayaUI/src/config.jl` (struct ~L26-29, `_build_config` ctor ~L107-109, `config_to_toml` ~L268-271)
- Test: `packages/HimalayaUI/test/test_config.jl`

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/test/test_config.jl`:

```julia
@testset "config round-trips beam center + pixel size" begin
    mktempdir() do dir
        toml = joinpath(dir, "experiment.toml")
        write(toml, """
        [experiment]
        name = "T/E"
        [beamline]
        energy_kev    = 12.0
        flight_path_m = 2.5
        beam_center_x = 420.791
        beam_center_y = 838.83
        pixel_size_um = 172.0
        [manifest]
        header_row = 0
        sample_id = 1
        name = 2
        display_name = 3
        filenames = 9
        notes_sample = 10
        notes_exposure = 11
        [layout]
        data_dir = "data"
        analysis_dir = "analysis/automatic_analysis"
        exposure_type = "simple"
        [files]
        integration = "{name}.dat"
        image = "{name}.tiff"
        """)
        cfg = HimalayaUI.load_config(toml)
        @test cfg.beam_center_x == 420.791
        @test cfg.beam_center_y == 838.83
        @test cfg.pixel_size_um == 172.0

        # Round-trip through config_to_toml preserves the populated fields
        # (this is the reingest persistence path — no mirror column).
        cfg2 = HimalayaUI._build_config(TOML.parse(HimalayaUI.config_to_toml(cfg)))
        @test cfg2.beam_center_x == 420.791
        @test cfg2.beam_center_y == 838.83
        @test cfg2.pixel_size_um == 172.0

        # Bare integer in TOML coerces to Float64.
        write(toml, replace(read(toml, String), "beam_center_x = 420.791" => "beam_center_x = 420"))
        @test HimalayaUI.load_config(toml).beam_center_x === 420.0
    end
end

@testset "config beam center absent round-trips to nothing" begin
    cfg = HimalayaUI.load_builtin_config("simple")
    @test cfg.beam_center_x === nothing
    @test cfg.beam_center_y === nothing
    @test cfg.pixel_size_um === nothing
    cfg2 = HimalayaUI._build_config(TOML.parse(HimalayaUI.config_to_toml(cfg)))
    @test cfg2.beam_center_x === nothing
    @test cfg2.pixel_size_um === nothing
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_config.jl`
Expected: FAIL — `type ExperimentConfig has no field beam_center_x` (or an `ExperimentConfig` constructor arity error).

- [ ] **Step 3: Add the struct fields**

In `packages/HimalayaUI/src/config.jl`, in `struct ExperimentConfig`, the `# [beamline]` block currently ends at `q_units ::String`. Add three fields immediately after it:

```julia
    # [beamline]
    energy_kev         ::Union{Float64,Nothing}
    flight_path_m      ::Union{Float64,Nothing}
    q_units            ::String
    beam_center_x      ::Union{Float64,Nothing}
    beam_center_y      ::Union{Float64,Nothing}
    pixel_size_um      ::Union{Float64,Nothing}
```

- [ ] **Step 4: Read them in `_build_config`**

In the `ExperimentConfig(...)` constructor call inside `_build_config`, immediately after the `get(bl, "q_units", "A-1"),` line, insert (order MUST match the struct):

```julia
        get(bl,  "q_units",       "A-1"),
        get(bl,  "beam_center_x", nothing),
        get(bl,  "beam_center_y", nothing),
        get(bl,  "pixel_size_um", nothing),
```

(The `Union{Float64,Nothing}` field auto-converts a bare `Int64` 420 to `420.0`.)

- [ ] **Step 5: Write them in `config_to_toml`**

In `config_to_toml`, after the `flight_path_m` omit-when-nothing line and before `beamline["q_units"] = cfg.q_units`, add:

```julia
    cfg.flight_path_m !== nothing && (beamline["flight_path_m"] = cfg.flight_path_m)
    cfg.beam_center_x !== nothing && (beamline["beam_center_x"] = cfg.beam_center_x)
    cfg.beam_center_y !== nothing && (beamline["beam_center_y"] = cfg.beam_center_y)
    cfg.pixel_size_um !== nothing && (beamline["pixel_size_um"] = cfg.pixel_size_um)
    beamline["q_units"] = cfg.q_units
```

- [ ] **Step 6: Run test to verify it passes**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_config.jl`
Expected: PASS (all testsets, including the two new ones).

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/src/config.jl packages/HimalayaUI/test/test_config.jl
git commit -m "feat(config): round-trip beam_center_x/y + pixel_size_um in [beamline]"
```

---

## Task 2: Experiment route — surface beam center + pixel size in JSON

**Files:**
- Modify: `packages/HimalayaUI/src/routes_experiments.jl` (`_q_units_from_config` L13-23, `_experiment_row_to_json` L25-29)
- Test: `packages/HimalayaUI/test/test_routes_experiments.jl`

Note: `_q_units_from_config` is **also called by `routes_samples.jl:53`** — do NOT remove it. Add `_beamline_from_config` and make `_q_units_from_config` delegate to it.

- [ ] **Step 1: Write the failing test**

In `packages/HimalayaUI/test/test_routes_experiments.jl`, the `@testset "experiments routes"` block already sets a malformed config and asserts `q_units == "A-1"`. Append a new testset at the end of the file:

```julia
@testset "experiment route surfaces beam center + pixel size" begin
    tmp = mktempdir()
    db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
    exp_id = HimalayaUI.create_experiment!(db; path=tmp, data_dir="data", analysis_dir="analysis")

    DBInterface.execute(db, "UPDATE experiments SET config = ? WHERE id = ?", [
        """
        [beamline]
        beam_center_x = 420.791
        beam_center_y = 838.83
        pixel_size_um = 172
        """, exp_id])

    HimalayaUI.with_test_server(db) do port, base
        body = JSON3.read(String(HTTP.get("$base/api/experiments/$exp_id").body))
        @test body.beam_center_x == 420.791
        @test body.beam_center_y == 838.83
        @test body.pixel_size_um == 172.0   # bare int coerced to Float64

        # Garbage numeric value must not 500 the route — falls back to null.
        DBInterface.execute(db, "UPDATE experiments SET config = ? WHERE id = ?",
            ["[beamline]\nbeam_center_x = \"oops\"\n", exp_id])
        r = HTTP.get("$base/api/experiments/$exp_id"; status_exception=false)
        @test r.status == 200
        @test JSON3.read(String(r.body)).beam_center_x === nothing
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (from `packages/HimalayaUI/`):
```bash
julia --project=. -e 'using Test, HimalayaUI, SQLite, HTTP, JSON3, DBInterface, Tables; include("test/test_http.jl"); include("test/test_routes_experiments.jl")'
```
Expected: FAIL — `type Object has no field beam_center_x` (the key is absent from the JSON).

- [ ] **Step 3: Add `_beamline_from_config` and delegate**

In `packages/HimalayaUI/src/routes_experiments.jl`, replace the `_q_units_from_config` function (L13-23) with:

```julia
"""
    _beamline_from_config(cfg_text) -> NamedTuple

Extract the render-only `[beamline]` fields from an experiment's stored config
TOML: `q_units` (String, default "A-1") and `beam_center_x`, `beam_center_y`,
`pixel_size_um` (each `Float64` or `nothing`). A bare integer in the TOML is
coerced to `Float64` to match the REAL-column behaviour of energy_kev.

Defensive by design: a non-string, empty, malformed, or non-numeric value all
fall back to the all-default tuple — one experiment's bad config must never 500
a list endpoint. The numeric coercion lives INSIDE the try so a value like
`beam_center_x = "oops"` throws into the fallback rather than out of the route.
"""
function _beamline_from_config(cfg_text)
    default = (q_units = "A-1", beam_center_x = nothing,
               beam_center_y = nothing, pixel_size_um = nothing)
    if cfg_text isa AbstractString && !isempty(cfg_text)
        try
            bl = get(TOML.parse(cfg_text), "beamline", Dict())
            num(k) = (v = get(bl, k, nothing); v === nothing ? nothing : Float64(v))
            return (q_units       = get(bl, "q_units", "A-1"),
                    beam_center_x = num("beam_center_x"),
                    beam_center_y = num("beam_center_y"),
                    pixel_size_um = num("pixel_size_um"))
        catch
            return default
        end
    end
    return default
end

# Back-compat shim: routes_samples.jl still calls this for its per-sample q_units.
_q_units_from_config(cfg_text)::String = _beamline_from_config(cfg_text).q_units
```

- [ ] **Step 4: Merge the fields in `_experiment_row_to_json`**

Replace `_experiment_row_to_json` (L25-29) with:

```julia
function _experiment_row_to_json(row::NamedTuple)
    d  = row_to_json(row)
    bl = _beamline_from_config(get(d, :config, nothing))
    d[:q_units]       = bl.q_units
    d[:beam_center_x] = bl.beam_center_x
    d[:beam_center_y] = bl.beam_center_y
    d[:pixel_size_um] = bl.pixel_size_um
    d
end
```

- [ ] **Step 5: Run test to verify it passes**

Run (from `packages/HimalayaUI/`):
```bash
julia --project=. -e 'using Test, HimalayaUI, SQLite, HTTP, JSON3, DBInterface, Tables; include("test/test_http.jl"); include("test/test_routes_experiments.jl")'
```
Expected: PASS (existing `q_units` assertions still green via the delegating shim; new testset green).

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/src/routes_experiments.jl packages/HimalayaUI/test/test_routes_experiments.jl
git commit -m "feat(routes): surface beam_center + pixel_size on experiment JSON"
```

---

## Task 3: Image route — emit raw detector dimensions as headers

**Files:**
- Modify: `packages/HimalayaUI/src/routes_exposures.jl` (image route, the `bytes = if is_thumb … else … end` block ~L58-74 and the `HTTP.Response` ~L80-84)
- Test: `packages/HimalayaUI/test/test_routes_image.jl`

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/test/test_routes_image.jl`:

```julia
@testset "full image emits raw X-Image-Width/Height; thumb does not" begin
    big_path = tempname() * ".tiff"
    save(big_path, Gray.(rand(Float32, 2048, 2048)))      # > 1536 → resized for display
    rect_path = tempname() * ".tiff"
    save(rect_path, Gray.(rand(Float32, 400, 600)))        # (rows=400, cols=600) → W=600, H=400

    db = SQLite.DB()
    HimalayaUI.create_schema!(db); HimalayaUI.migrate_schema!(db)
    exp  = HimalayaUI.create_experiment!(db; path="/tmp", data_dir="/tmp", analysis_dir="/tmp")
    samp = HimalayaUI.create_sample!(db; experiment_id=exp)
    big  = HimalayaUI.create_exposure!(db; sample_id=samp, image_path=big_path)
    rect = HimalayaUI.create_exposure!(db; sample_id=samp, image_path=rect_path)

    HimalayaUI.with_test_server(db) do port, base
        rb = HTTP.get("$base/api/exposures/$big/image")
        h  = Dict(rb.headers)
        @test h["X-Image-Width"]  == "2048"   # RAW, not the <=1536 displayed size
        @test h["X-Image-Height"] == "2048"

        rr = HTTP.get("$base/api/exposures/$rect/image")
        hr = Dict(rr.headers)
        @test hr["X-Image-Width"]  == "600"   # cols
        @test hr["X-Image-Height"] == "400"   # rows

        rt = HTTP.get("$base/api/exposures/$big/image?thumb=1")
        @test !haskey(Dict(rt.headers), "X-Image-Width")
    end
    rm(big_path; force=true); rm(rect_path; force=true)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (from `packages/HimalayaUI/`):
```bash
julia --project=. -e 'using Test, HimalayaUI, SQLite, HTTP, JSON3, DBInterface, FileIO, ImageCore, TiffImages; include("test/test_http.jl"); include("test/test_routes_image.jl")'
```
Expected: FAIL — `KeyError: key "X-Image-Width" not found`.

- [ ] **Step 3: Capture raw dims + emit headers**

In `packages/HimalayaUI/src/routes_exposures.jl`, the image route. Before the `bytes = if is_thumb` line, declare the raw-dim holders, capture them in the full branch, and build the response headers conditionally. Replace from `bytes = if is_thumb` through the `HTTP.Response(...)` call with:

```julia
        raw_w = 0; raw_h = 0
        bytes = if is_thumb
            ensure_thumb_cached(db, id, path; token = vtoken)
        else
            # Full path: percentile-clip at full res, capture RAW dims (the beam
            # center's coordinate space) BEFORE the 1536 display cap, then resize.
            img = load_and_lognormalize(path)
            raw_h, raw_w = size(img)          # Julia size = (rows, cols) = (h, w)
            img = resize_to_fit(img, 1536)
            encode_png(img)
        end

        hdrs = ["Content-Type"    => "image/png",
                "Cache-Control"   => "private, max-age=31536000, immutable",
                "X-Image-Version" => vtoken]
        if !is_thumb
            # Raw detector pixel dims for the q-ring overlay calibration. Full
            # branch only — the thumb never carries rings. Frontend tolerates
            # their absence (→ centered fallback).
            push!(hdrs, "X-Image-Width"  => string(raw_w))
            push!(hdrs, "X-Image-Height" => string(raw_h))
        end
        HTTP.Response(200, hdrs, bytes)
```

(Keep the existing comments above the `if is_thumb` branches if you prefer; the key changes are the `raw_w/raw_h` capture and the conditional `hdrs`.)

- [ ] **Step 4: Run test to verify it passes**

Run (from `packages/HimalayaUI/`):
```bash
julia --project=. -e 'using Test, HimalayaUI, SQLite, HTTP, JSON3, DBInterface, FileIO, ImageCore, TiffImages; include("test/test_http.jl"); include("test/test_routes_image.jl")'
```
Expected: PASS (new testset + the existing 1536-cap and 404 testsets).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_exposures.jl packages/HimalayaUI/test/test_routes_image.jl
git commit -m "feat(routes): emit raw X-Image-Width/Height on full detector image"
```

---

## Task 4: Frontend types + `buildDetectorCalibration` helper

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/api.ts` (`interface Experiment`, after `q_units` ~L16)
- Modify: `packages/HimalayaUI/frontend/src/print/pages/focusAdapters.ts` (new export; imports near top)
- Test: `packages/HimalayaUI/frontend/test/print-pages/focusAdapters.test.ts`

- [ ] **Step 1: Write the failing test**

Add the two imports to the **top** of `packages/HimalayaUI/frontend/test/print-pages/focusAdapters.test.ts` (alongside the existing imports), then append the `describe` block at the end of the file:

```ts
// add to the existing import group at the top:
import { buildDetectorCalibration } from "../../src/print/pages/focusAdapters";
import type { Experiment } from "../../src/api";
```

```ts
// append at the end of the file:
function makeExp(over: Partial<Experiment> = {}): Experiment {
  return {
    id: 1, name: "E", path: "/p", data_dir: "d", analysis_dir: "a",
    manifest_path: null, created_at: "", q_units: "A-1",
    beam_center_x: 420.791, beam_center_y: 838.83, pixel_size_um: 172,
    energy_kev: 12, flight_path_m: 2.5, ...over,
  } as Experiment;
}

describe("buildDetectorCalibration", () => {
  it("builds a full calibration with µm→mm and m→mm conversions", () => {
    const cal = buildDetectorCalibration(makeExp(), { w: 1475, h: 1679 });
    expect(cal).toEqual({
      beamCenterPx: { x: 420.791, y: 838.83 },
      imageSizePx: { w: 1475, h: 1679 },
      sampleDistanceMm: 2500,
      pixelSizeMm: 0.172,
      energyKeV: 12,
    });
  });

  it("returns null when rawSize is missing", () => {
    expect(buildDetectorCalibration(makeExp(), null)).toBeNull();
  });

  it("returns null when experiment is undefined", () => {
    expect(buildDetectorCalibration(undefined, { w: 100, h: 100 })).toBeNull();
  });

  it.each(["beam_center_x", "beam_center_y", "pixel_size_um", "energy_kev", "flight_path_m"] as const)(
    "returns null when %s is null", (field) => {
      expect(buildDetectorCalibration(makeExp({ [field]: null }), { w: 100, h: 100 })).toBeNull();
    });

  it("returns null for non-positive or non-finite dims", () => {
    expect(buildDetectorCalibration(makeExp(), { w: 0, h: 100 })).toBeNull();
    expect(buildDetectorCalibration(makeExp(), { w: NaN, h: 100 })).toBeNull();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- focusAdapters`
Expected: FAIL — `buildDetectorCalibration is not a function` (and a tsc error on the new `Experiment` fields).

- [ ] **Step 3: Add the `Experiment` fields**

In `packages/HimalayaUI/frontend/src/api.ts`, in `interface Experiment`, replace the `q_units` line with:

```ts
  q_units: string | null;
  beam_center_x: number | null;
  beam_center_y: number | null;
  pixel_size_um: number | null;
  energy_kev: number | null;
  flight_path_m: number | null;
```

- [ ] **Step 4: Add the helper**

In `packages/HimalayaUI/frontend/src/print/pages/focusAdapters.ts`, add imports near the existing ones at the top:

```ts
import type { Experiment } from "../../api";
import type { DetectorCalibration } from "../detector/detectorGeometry";
```

Then add this exported function (place it next to `toDetectorRings`, the other ring-feeding adapter):

```ts
/**
 * Experiment beamline params + the RAW detector pixel size (from the image
 * route's X-Image-Width/Height headers) → a DetectorCalibration for the
 * geometry engine. Returns null unless ALL ingredients are present and finite
 * (a 0/NaN dim or null field would yield NaN/Infinity radii); null → the
 * DetectorPanel centered fallback. Pure + unit-tested; all the arithmetic and
 * guarding lives here, not in the component.
 */
export function buildDetectorCalibration(
  experiment: Experiment | undefined,
  rawSize: { w: number; h: number } | null,
): DetectorCalibration | null {
  if (!experiment || !rawSize) return null;
  const { beam_center_x, beam_center_y, pixel_size_um, energy_kev, flight_path_m } = experiment;
  const ok = (n: number | null): n is number => n !== null && Number.isFinite(n);
  if (!ok(beam_center_x) || !ok(beam_center_y) || !ok(pixel_size_um) ||
      !ok(energy_kev) || !ok(flight_path_m)) return null;
  if (!(rawSize.w > 0) || !(rawSize.h > 0)) return null;
  return {
    beamCenterPx: { x: beam_center_x, y: beam_center_y },
    imageSizePx: { w: rawSize.w, h: rawSize.h },
    sampleDistanceMm: flight_path_m * 1000,
    pixelSizeMm: pixel_size_um / 1000,
    energyKeV: energy_kev,
  };
}
```

- [ ] **Step 5: Run test to verify it passes**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- focusAdapters`
Expected: PASS (all cases).

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/api.ts packages/HimalayaUI/frontend/src/print/pages/focusAdapters.ts packages/HimalayaUI/frontend/test/print-pages/focusAdapters.test.ts
git commit -m "feat(focus): add buildDetectorCalibration + Experiment beamline fields"
```

---

## Task 5: `DetectorImage` — surface raw size + orientation via callbacks

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/print/detector/DetectorImage.tsx` (`Props` L5-15, `evaluateOrient` L28-38, `renderImage` L40-69, the component body for refs)
- Test: `packages/HimalayaUI/frontend/test/print-detector/DetectorImage.test.tsx`

- [ ] **Step 1: Write the failing test (and fix the existing mock)**

In `packages/HimalayaUI/frontend/test/print-detector/DetectorImage.test.tsx`, update the `beforeEach` `fetch` mock to include a `headers.get` stub:

```ts
  global.fetch = vi.fn().mockResolvedValue({
    ok: true,
    headers: { get: (k: string) => (k === "X-Image-Width" ? "2048" : k === "X-Image-Height" ? "1024" : null) },
    blob: () => Promise.resolve(new Blob(
      [Uint8Array.from(atob(TINY_PNG), (c) => c.charCodeAt(0))], { type: "image/png" })),
  } as unknown as Response);
```

Then append two tests:

```ts
test("reports raw image size from X-Image-Width/Height headers", async () => {
  const onRawSize = vi.fn();
  render(<DetectorImage src="/x.png" size="full" onRawSize={onRawSize} />);
  await waitFor(() => expect(onRawSize).toHaveBeenCalledWith(2048, 1024));
});

test("a headerless response does not call onRawSize and does not throw", async () => {
  global.fetch = vi.fn().mockResolvedValue({
    ok: true,
    blob: () => Promise.resolve(new Blob([new Uint8Array(0)], { type: "image/png" })),
  } as unknown as Response);
  const onRawSize = vi.fn();
  render(<DetectorImage src="/x.png" size="full" onRawSize={onRawSize} />);
  await waitFor(() => expect(screen.getByRole("img", { hidden: true })).toBeInTheDocument());
  expect(onRawSize).not.toHaveBeenCalled();
});
```

- [ ] **Step 2: Run test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- DetectorImage`
Expected: FAIL — `onRawSize` not called (prop doesn't exist / headers not read).

- [ ] **Step 3: Add props + callback refs + read headers + orient transition**

In `packages/HimalayaUI/frontend/src/print/detector/DetectorImage.tsx`:

(a) Extend `Props`:

```ts
interface Props {
  src: string | null;
  size: "thumb" | "full";
  lutVariant?: DetectorLutVariant;
  className?: string;
  /** Raw detector pixel size from X-Image-Width/Height headers (notify-only). */
  onRawSize?: (w: number, h: number) => void;
  /** Display orientation, fired only on a true transition (notify-only). */
  onOrient?: (o: "portrait" | "landscape") => void;
}
```

(b) Destructure them and add stable refs at the top of the component (so they are NOT render-input deps of `renderImage`/`evaluateOrient`):

```ts
export function DetectorImage({ src, size, lutVariant = "neutral", className, onRawSize, onOrient }: Props): JSX.Element {
  const wrapperRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const onRawSizeRef = useRef(onRawSize); onRawSizeRef.current = onRawSize;
  const onOrientRef = useRef(onOrient); onOrientRef.current = onOrient;
  const rawSizeRef = useRef<{ w: number; h: number } | null>(null);
  const prevOrientRef = useRef<"portrait" | "landscape">("portrait");
```

(c) Rewrite `evaluateOrient` to compute one `next` layout, set it, and fire `onOrient` only on a true transition (it previously had two `setLayout` sites):

```ts
  const evaluateOrient = useCallback(() => {
    const wrapper = wrapperRef.current, canvas = canvasRef.current;
    if (!wrapper || !canvas || !canvas.width || !canvas.height) return;
    let next: Layout;
    if (size === "thumb") {
      next = { orient: "portrait", caps: null };
    } else {
      const cw = wrapper.clientWidth, ch = wrapper.clientHeight;
      if (cw === 0 || ch === 0) return;
      next = decideOrient({
        containerW: cw, containerH: ch, imageW: canvas.width, imageH: canvas.height,
        viewportW: typeof window !== "undefined" ? window.innerWidth : 0,
      });
    }
    setLayout(next);
    if (prevOrientRef.current !== next.orient) {
      prevOrientRef.current = next.orient;
      onOrientRef.current?.(next.orient);
    }
  }, [size]);
```

(d) In `renderImage`, right after `if (!res.ok) return;`, read the headers and report a changed raw size:

```ts
    const res = await fetch(src);
    if (!res.ok) return;
    const rw = Number(res.headers?.get?.("X-Image-Width"));
    const rh = Number(res.headers?.get?.("X-Image-Height"));
    if (Number.isFinite(rw) && Number.isFinite(rh) && rw > 0 && rh > 0 &&
        (rawSizeRef.current?.w !== rw || rawSizeRef.current?.h !== rh)) {
      rawSizeRef.current = { w: rw, h: rh };
      onRawSizeRef.current?.(rw, rh);
    }
```

Leave `renderImage`'s dep array as `[src, hasIntersected, lutVariant, evaluateOrient]` — do NOT add the callbacks (they are accessed via refs).

- [ ] **Step 4: Run test to verify it passes**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- DetectorImage`
Expected: PASS — including the existing canvas-render / LUT / object-fit / orient tests (the `transform === ""` portrait test still holds; `headers` stub is benign).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/detector/DetectorImage.tsx packages/HimalayaUI/frontend/test/print-detector/DetectorImage.test.tsx
git commit -m "feat(detector): DetectorImage reports raw size + orient via callbacks"
```

---

## Task 6: `DetectorPanel` — forward callbacks + pass orient to rings

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/print/components/DetectorPanel.tsx` (`DetectorPanelProps` L13-52, body destructure L54-68, `DetectorImage` call L96-101, `DetectorRings` call L102-110)
- Test: `packages/HimalayaUI/frontend/test/print-components/DetectorPanel.test.tsx`

- [ ] **Step 1: Write the failing test**

Append to `packages/HimalayaUI/frontend/test/print-components/DetectorPanel.test.tsx`:

```ts
test("passes orient through to DetectorRings", async () => {
  render(<DetectorPanel src="/x.png" rings={[0.1, 0.2]} orient="landscape" />);
  const rings = await screen.findByTestId("detector-rings");
  expect(rings.getAttribute("data-orient")).toBe("landscape");
});

test("forwards onRawSize/onOrient to DetectorImage (props accepted, no crash)", () => {
  const onRawSize = vi.fn(); const onOrient = vi.fn();
  expect(() =>
    render(<DetectorPanel src="/x.png" onRawSize={onRawSize} onOrient={onOrient} />),
  ).not.toThrow();
});
```

(If `DetectorPanel.test.tsx` does not already mock `fetch`/`createImageBitmap`, copy the `beforeEach` mock block from `DetectorImage.test.tsx` — including the `headers.get` stub — so the canvas renderer resolves.)

- [ ] **Step 2: Run test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- DetectorPanel`
Expected: FAIL — tsc rejects the unknown `orient`/`onRawSize`/`onOrient` props (and `data-orient` would be `"portrait"`).

- [ ] **Step 3: Add the props**

In `packages/HimalayaUI/frontend/src/print/components/DetectorPanel.tsx`, add to `DetectorPanelProps` (after `imageAspect?`):

```ts
  /** Display orientation, forwarded to DetectorRings so the overlay rotates in
   *  lockstep with a landscape-rotated canvas. */
  orient?: "portrait" | "landscape";
  /** Raw detector pixel size from the image headers (notify-only). */
  onRawSize?: (w: number, h: number) => void;
  /** Display-orientation transition (notify-only). */
  onOrient?: (o: "portrait" | "landscape") => void;
```

Add `orient, onRawSize, onOrient` to the destructured params in the function signature.

- [ ] **Step 4: Forward them**

In the `DetectorImage` JSX, add the callback forwards (conditional spread for `exactOptionalPropertyTypes`):

```tsx
        <DetectorImage
          src={src}
          size="full"
          className="h-full w-full"
          {...(lutVariant !== undefined ? { lutVariant } : {})}
          {...(onRawSize !== undefined ? { onRawSize } : {})}
          {...(onOrient !== undefined ? { onOrient } : {})}
        />
```

In the `DetectorRings` JSX, add the `orient` pass-through:

```tsx
          <DetectorRings
            beamCenter={center}
            rings={placed.rings}
            imageAspect={aspect}
            {...(orient !== undefined ? { orient } : {})}
            {...(hoveredQ !== undefined ? { hoveredQ } : {})}
            {...(onHoverQ !== undefined ? { onHoverQ } : {})}
          />
```

- [ ] **Step 5: Run test to verify it passes**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- DetectorPanel`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/components/DetectorPanel.tsx packages/HimalayaUI/frontend/test/print-components/DetectorPanel.test.tsx
git commit -m "feat(detector): DetectorPanel forwards orient + size callbacks"
```

---

## Task 7: `FocusPage` — wire calibration, imageAspect, orient

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/print/pages/FocusPage.tsx` (imports; component state near `experimentQ` L163; the main `DetectorPanel` JSX L785-806)
- Verification: tsc (`npm run build` in Task 9) + existing FocusPage suite; behavioral confirmation in Task 9 live verify.

This task has no jsdom unit test: the calibration only materializes after the canvas `fetch` resolves the header size, which is integration-level. The pure logic is already unit-tested (Task 4); correctness of the wiring is confirmed by the live render step (Task 9). Do not write a fake/over-mocked test.

- [ ] **Step 1: Add the import**

In `FocusPage.tsx`, add `buildDetectorCalibration` to the existing import from `./focusAdapters` (or its own import line):

```ts
import { buildDetectorCalibration } from "./focusAdapters";
```

- [ ] **Step 2: Add state + derived calibration**

Near the other FocusPage state (after `const experimentQ = useExperiment(...)`, ~L163), add:

```ts
  const [rawSize, setRawSize] = useState<{ w: number; h: number } | null>(null);
  const [detOrient, setDetOrient] = useState<"portrait" | "landscape">("portrait");
  const handleRawSize = useCallback(
    (w: number, h: number) => setRawSize((p) => (p?.w === w && p?.h === h ? p : { w, h })),
    [],
  );
  const calibration = buildDetectorCalibration(experimentQ.data, rawSize);
  const imageAspect = rawSize ? rawSize.w / rawSize.h : undefined;
```

(Confirm `useState`/`useCallback` are imported from `react` at the top — add them if missing.)

- [ ] **Step 3: Pass to the main DetectorPanel**

In the main `DetectorPanel` JSX (the one with `src={detectorSrc}`, ~L785), add the calibration + orientation props (conditional spread for the optional ones):

```tsx
              <DetectorPanel
                src={detectorSrc}
                rings={rings}
                ringPhases={ringPhases}
                orient={detOrient}
                onRawSize={handleRawSize}
                onOrient={setDetOrient}
                {...(calibration ? { calibration } : {})}
                {...(imageAspect !== undefined ? { imageAspect } : {})}
                {...(hoveredQ !== undefined ? { hoveredQ } : {})}
                onHoverQ={setHoveredQ}
                tools={ /* unchanged */ }
              />
```

(Leave the `tools={…}` block exactly as it is — only add the new props.)

- [ ] **Step 4: Typecheck + existing suite**

Run (from `packages/HimalayaUI/frontend/`):
```bash
npm test -- FocusPage
npx tsc --noEmit -p tsconfig.json
```
Expected: existing FocusPage tests PASS; tsc reports no errors.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/pages/FocusPage.tsx
git commit -m "feat(focus): wire beam-center calibration + imageAspect into DetectorPanel"
```

---

## Task 8: Template + docs

**Files:**
- Modify: `packages/HimalayaUI/configs/simple.toml` (`[beamline]` block)
- Modify: `docs/experiment-config.md` (`[beamline]` table ~L101-110)

- [ ] **Step 1: Update the template**

In `packages/HimalayaUI/configs/simple.toml`, replace the `[beamline]` block's commented lines with:

```toml
[beamline]
# Optional. Stored as NULL when omitted (preferred over a 0.0 sentinel,
# which would misrepresent "unknown" as a real measurement of zero).
# Uncomment and set when you know them:
# energy_kev    = 12.0    # X-ray photon energy in keV
# flight_path_m =  1.5    # sample-to-detector distance in metres
# q_units       = "A-1"  # axis label units (ASCII; UI prettifies to Å⁻¹). Alternatives: "nm-1", "1/A"
#
# Beam center (detector pixels, origin BOTTOM-LEFT, y-up) + pixel pitch enable
# physically-correct q-rings on the Focus page. All five beamline values
# (energy_kev, flight_path_m, beam_center_x, beam_center_y, pixel_size_um) must
# be set; omit any and the rings fall back to a centered, decorative overlay.
# beam_center_x = 420.791
# beam_center_y = 838.83
# pixel_size_um = 172.0   # detector pixel pitch in microns (e.g. Pilatus 172)
```

- [ ] **Step 2: Update the docs**

In `docs/experiment-config.md`, replace the `[beamline]` field table (the rows for `energy_kev`/`flight_path_m`) with:

```markdown
| Field | Purpose |
|-------|---------|
| `energy_kev` | Beam energy in keV. Used for q-ring radii (q-calibration). |
| `flight_path_m` | Sample-to-detector distance in metres. Used for q-ring radii. |
| `beam_center_x` | Beam center column, in **detector pixels, origin bottom-left** (y-up). |
| `beam_center_y` | Beam center row, in **detector pixels, origin bottom-left** (y-up). |
| `pixel_size_um` | Detector pixel pitch in microns (e.g. Pilatus 172, Eiger 75). |

The Focus page draws physically-correct, beam-centered q-rings only when **all
five** beamline values are present (`energy_kev`, `flight_path_m`,
`beam_center_x`, `beam_center_y`, `pixel_size_um`). Omit any and the rings fall
back to a centered, decorative overlay — no error. Beam center and pixel size
are render-only: they live in the config blob, are surfaced on the experiment
API, and (unlike `energy_kev`/`flight_path_m`) are not mirrored to queryable DB
columns.
```

- [ ] **Step 3: Verify the template still loads**

Run: `julia --project=packages/HimalayaUI packages/HimalayaUI/test/test_config.jl`
Expected: PASS (the `load_builtin_config simple` testset proves the edited template parses; the commented beam-center lines stay `nothing`).

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/configs/simple.toml docs/experiment-config.md
git commit -m "docs: document [beamline] beam center + pixel size for q-rings"
```

---

## Task 9: Full gates + live render verification

**Files:** none (verification only).

- [ ] **Step 1: Full frontend gates**

Run (from `packages/HimalayaUI/frontend/`):
```bash
npm test
npm run build        # tsc --noEmit + vite build
npm run lint:design  # design-system guard (new props are data/geometry, not appearance)
```
Expected: all green. The new `onRawSize`/`onOrient`/`orient`/`calibration`/`imageAspect` props are data callbacks + geometry and must not trip `lint:design`.

- [ ] **Step 2: Full Julia suite (capture once)**

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' > /tmp/jl-test.out 2>&1
grep -E "Test Summary|did not pass|fail" /tmp/jl-test.out ; tail -50 /tmp/jl-test.out
```
Expected: no failures. (Confirms `routes_samples` q_units still works via the delegating shim, and the read-only-dir + reingest regressions stay green.)

- [ ] **Step 3: Live render verification**

Serve a real experiment that has a known beam center in its `experiment.toml`, open a sample's Focus page, and confirm with the browser (see `docs/frontend-dev-loop.md` + the live-audit harness):
- The q-rings are centered on the beamstop (not the image middle).
- A known-q ring lands on the corresponding diffraction radius.
- At viewport ≥ 1400px with a wide detector column, there is no flicker-rotate on first paint.

**Disambiguate a vertical-mirror symptom** if rings look flipped: (a) wrong beam convention → flip is a one-line change in `buildRingPlacements` (`y = 1 - …`); (b) the beamline writes bottom-up TIFFs so the *image* is flipped relative to the beam space. Confirm which before "fixing".

- [ ] **Step 4: Final commit (if any verification tweaks)**

```bash
git add -A
git commit -m "test: full-suite + live-render verification for beam-center calibration"
```

(Only if a convention flip or test tweak was needed; otherwise nothing to commit.)

---

## Self-Review notes (author)

- **Spec coverage:** §1→T1, §2→T2, §2b→T3, §3→T4, §4→T4, §5→T5, §6(panel)→T6, §6(page)→T7, §7→T8, Testing+live→T1–T3/T4/T9. All spec sections mapped.
- **Type consistency:** `buildDetectorCalibration(experiment, rawSize)` signature identical in T4 (def) and T7 (call); `DetectorCalibration` field names (`beamCenterPx`, `imageSizePx`, `sampleDistanceMm`, `pixelSizeMm`, `energyKeV`) match `detectorGeometry.ts`; `onRawSize(w,h)`/`onOrient(o)` identical across T5/T6/T7; `X-Image-Width`/`X-Image-Height` header names identical in T3 (emit) and T5 (read).
- **Known judgment call:** T7 has no jsdom unit test by design (the calibration is gated on an async canvas fetch + header read; the pure logic is covered in T4, behavior in T9 live verify). Flagged explicitly rather than writing an over-mocked test.
- **No-removal note:** `_q_units_from_config` is retained (delegating shim) because `routes_samples.jl:53` still calls it.
