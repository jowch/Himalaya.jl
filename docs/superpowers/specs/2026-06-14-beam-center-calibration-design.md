# Beam-center calibration for detector q-rings

**Status:** design approved, awaiting spec review
**Date:** 2026-06-14

## Problem

The Focus page overlays phase-coloured Debye–Scherrer (q-)rings on the
detector image. The rendering engine (`detectorGeometry.ts` +
`DetectorRings`) is already complete and can place rings either:

- **fully calibrated** — beam center + physically-correct radii, given a
  `DetectorCalibration { beamCenterPx, imageSizePx, sampleDistanceMm,
  pixelSizeMm, energyKeV }`, or
- **centered fallback** — centered at the image middle (0.5, 0.5) with
  radii merely spread across the q-range (decorative).

Today `FocusPage` passes **neither** a `calibration` nor a `beamCenter`,
so the rings always land in the centered/decorative fallback. A ring
labelled `q = 0.1` does not sit on the actual `q = 0.1` diffraction
radius, and the ring origin ignores the true beam position.

There is no way to tell Himalaya where the beam center is, nor the
detector pixel size, so the existing full-calibration path is unreachable.

## Goal

Let a scientist declare the beam center and pixel size in
`experiment.toml`, thread those (with the already-present beam energy and
flight path) into the existing calibration engine, and render
**physically-correct, beam-centered** q-rings on the Focus page.

When the calibration is incomplete, behaviour is unchanged (centered
fallback) — no regression for existing experiments.

## Decisions (settled during brainstorming)

- **Full geometric calibration**, not center-only. Rings must be both
  centered on the true beam *and* truly sized.
- **Beam center convention:** detector pixels, origin **bottom-left**
  (physics y-up). This matches the assumption already baked into
  `buildRingPlacements` (it flips y: `y = 1 - beamCenterPx.y /
  imageSizePx.h`). Example value: `(420.791, 838.83)`.
  - Marked "probably bottom-left" by the user → **verify by rendering**
    against a real detector image; if rings land mirrored vertically the
    flip is a one-line change.
- **Pixel size unit:** microns (`pixel_size_um`) — how detector specs are
  quoted (Pilatus 172 µm, Eiger 75 µm). Converted to mm internally.
- **No DB migration / no mirror columns.** Beam center and pixel size
  live only in the `experiments.config` TOML blob and are parsed on read,
  following the **`q_units` precedent** (`_q_units_from_config`), not the
  `energy_kev`/`flight_path_m` mirror-column precedent. They are never
  queried, only rendered.
- **Raw detector dimensions come from the image-route response headers**,
  not config and not the displayed bitmap. Two facts force this:
  - `DetectorImage` does **not** use an `<img>` element — it `fetch`es the
    PNG, decodes it with `createImageBitmap`, applies a LUT on an
    offscreen canvas, and paints a `<canvas>` (DetectorImage.tsx). There
    is no `naturalWidth`/`onLoad`; the decoded bitmap's `width/height` is
    the only client-side size.
  - The `size="full"` route **downscales** the rendered PNG to ≤ 1536px
    max-side (`resize_to_fit(img, 1536)`, issue #260). So the displayed
    bitmap is **not** raw detector resolution for any detector larger than
    1536px. The beam center is in *raw* detector pixels, and the q→radius
    conversion also divides by the raw width, so both need the **raw**
    dimensions — using the resized bitmap size would inflate the center
    offset and every radius by `1/scale`.
  - `resize_to_fit` scales uniformly, so the displayed bitmap's *aspect*
    is still correct (used for the frame box and the ring y-correction);
    only its absolute size is wrong.
  - Fix: the image route emits `X-Image-Width` / `X-Image-Height` headers
    carrying the **pre-resize** raster size (`size(img)` before
    `resize_to_fit`). `DetectorImage` reads them off the `fetch` response
    it already makes and surfaces them via callback. Normalized fractions
    (beam center, radii) are scale-invariant, so computing them against
    the raw size renders correctly on the downscaled display.
- **Landscape auto-rotation must be handled.** `DetectorImage` rotates the
  canvas 90° for wide containers (`decideOrient`, gated to viewport ≥
  1400px when `containerAspect > imageAspect × 1.25`). `DetectorRings`
  already accepts an `orient` prop and rotates the SVG in lockstep, but
  `DetectorPanel` never passes it — invisible today only because centered
  concentric circles are rotation-invariant. An **off-center** beam would
  desync rings from the image. Two-part fix:
  - Pass the real `imageAspect` (see below) so the frame box matches the
    image aspect; then `containerAspect ≈ imageAspect`, the rotation
    predicate is false, and rotation does not trigger for the detector
    panel.
  - Defensively surface `orient` from `DetectorImage` and pass it to
    `DetectorRings` so the overlay stays registered even if a future
    layout does trigger rotation.

## Config schema — `[beamline]`

Three new **optional** fields:

```toml
[beamline]
energy_kev     = 12.0       # existing — used for ring radii (q-calibration)
flight_path_m  = 2.5        # existing — sample-to-detector distance
beam_center_x  = 420.791    # NEW: detector pixels, origin bottom-left
beam_center_y  = 838.83     # NEW: detector pixels, origin bottom-left (y-up)
pixel_size_um  = 172.0      # NEW: detector pixel pitch in microns
```

The four ingredients required for full calibration are `beam_center_x`,
`beam_center_y`, `pixel_size_um`, and the pre-existing `energy_kev` +
`flight_path_m`. If **any** is absent, the Focus page renders the
centered fallback exactly as it does today.

## Components and changes

### 1. `packages/HimalayaUI/src/config.jl`

- Add to `struct ExperimentConfig`: `beam_center_x`, `beam_center_y`,
  `pixel_size_um`, each `::Union{Float64,Nothing}`.
- `_build_config`: read them from the `[beamline]` table with
  `get(bl, "beam_center_x", nothing)` etc.
- `config_to_toml`: write each into the `beamline` dict only when
  `!== nothing` (omit-when-unset, exactly like `energy_kev` /
  `flight_path_m`), so a round-trip preserves `nothing`. **This is the
  persistence path, not just a test concern:** `init` and
  `_reingest_inner!` (cli.jl) store `config_to_toml(cfg)` — the
  *re-serialized* blob, never the raw file text. Because there is no
  mirror column, if `config_to_toml` fails to write a populated field,
  reingest **silently strips the beam center with no trace**. The
  config round-trip test (see Testing) must assert the *populated* case
  survives, not only absent→`nothing`.

### 2. Backend route — `packages/HimalayaUI/src/routes_experiments.jl`

- Generalize `_q_units_from_config(cfg_text) -> String` into a single
  blob extractor that returns `q_units`, `beam_center_x`,
  `beam_center_y`, and `pixel_size_um` from **one** `TOML.parse`.
- `_experiment_row_to_json` merges those four into the experiment JSON.
  `energy_kev` / `flight_path_m` already ride along via `SELECT *`.
- Keep the defensive contract: a non-string, empty, or malformed config
  must never 500 a list endpoint — missing/garbage values resolve to
  `q_units = "A-1"` and `nothing` for the numeric fields. **The numeric
  coercion must live inside the single `try`/`catch`:** coerce each field
  to `Float64` (a bare `beam_center_x = 420` parses as `Int64`, and the
  JSON should be a float to match the `REAL`-column behaviour of
  `energy_kev`/`flight_path_m`), with `x === nothing ? nothing :
  Float64(x)`. A non-numeric garbage value (`beam_center_x = "oops"`)
  then throws *inside* the try and falls back to the all-default tuple,
  rather than 500-ing the list endpoint. Placing the coercion outside the
  try would break the contract.

### 2b. Image route — `packages/HimalayaUI/src/routes_exposures.jl`

The full-image branch loads the raster (`img = load_and_lognormalize(path)`)
and then downscales it (`resize_to_fit(img, 1536)`). Capture the
**pre-resize** dimensions (`h, w = size(img)` — Julia `size` is
`(rows, cols) = (height, width)`) and emit them as response headers
`X-Image-Width => w`, `X-Image-Height => h` alongside the existing
`X-Image-Version`. These are the raw detector pixel dimensions the beam
center is expressed in.

- **Full branch only.** `bytes` is assigned in two branches (thumb vs
  full); the raw dims exist only in the full branch. The thumb branch is
  unaffected — it never carries a ring overlay. The frontend must tolerate
  the headers' absence (→ centered fallback), never assume them.
- `load_and_lognormalize` applies **no transpose/rotation** (verified), so
  `size(img)` is the native TIFF raster orientation — consistent with the
  displayed canvas.
- The response is `Cache-Control: immutable` + `?v=<token>` URL-keyed, so
  the new headers ride the cached PNG on repeat loads (they are stable per
  raster — fine). The frontend must not assume a fresh network round-trip
  delivers them each time.

### 3. Frontend types — `packages/HimalayaUI/frontend/src/api.ts`

Add to `interface Experiment` (all `number | null`):
`beam_center_x`, `beam_center_y`, `pixel_size_um`, `energy_kev`,
`flight_path_m`. (The last two already flow in the JSON but are
currently untyped.)

### 4. Calibration assembly — `focusAdapters.ts` (new pure helper)

```ts
buildDetectorCalibration(
  experiment: Experiment | undefined,
  rawSize: { w: number; h: number } | null,
): DetectorCalibration | null
```

- `rawSize` is the **raw** detector pixel size from the image headers
  (§2b/§5), used as `imageSizePx`. Returns `null` if `experiment`,
  `rawSize`, or **any** of the four ingredients is missing/null.
- **Finite/positive guard.** Because `qToImageRadius` and
  `buildRingPlacements` divide by `imageSizePx.w`/`.h`, a header that
  parses to `0`/`NaN` (`Number("")` → `NaN`) would yield `NaN`/`Infinity`
  radii and silently render broken rings. Reject with `null` unless
  `Number.isFinite` holds for all five numeric inputs **and**
  `rawSize.w > 0 && rawSize.h > 0`. This is part of the helper contract,
  tested.
- Unit conversions: `sampleDistanceMm = flight_path_m * 1000`,
  `pixelSizeMm = pixel_size_um / 1000`.
- Maps `beam_center_x/y` → `beamCenterPx`, `rawSize` → `imageSizePx`,
  `energy_kev` → `energyKeV`.

This is the unit-tested seam; all arithmetic and null-guarding lives here,
not in the component. Co-locate it with `toDetectorRings` (the other
ring-feeding adapter) in `focusAdapters.ts`.

### 5. Image dimensions + orientation — `DetectorImage` / `DetectorPanel`

- `DetectorImage` (canvas renderer, no `<img>`): in `renderImage`, after
  the `fetch`, read the **raw** pre-resize size via
  `res.headers?.get?.("X-Image-Width" / "X-Image-Height")` — **defensive
  optional chaining is required**: the existing `DetectorImage` test mocks
  return `{ ok, blob }` with no `headers`, and a headerless/older cached
  response must degrade to the centered fallback, not throw. Report it via
  a new optional `onRawSize?: (w, h) => void`. Also report orientation via
  a new optional `onOrient?: (o: "portrait" | "landscape") => void`.
- **Callback-loop hazards (load-bearing):**
  - `onRawSize`/`onOrient` must be **excluded** from `renderImage`'s and
    `evaluateOrient`'s `useCallback`/effect dep arrays — they are
    notify-only outputs, not render inputs. Including them (with inline
    arrows from the parent) re-fires the fetch/orient effects every render
    → refetch/setState loop.
  - `onOrient` must fire only on a true **transition** (diff against a
    `useRef` of the previous orient), not on every `evaluateOrient` —
    `evaluateOrient` runs on every ResizeObserver tick and calls
    `setLayout` with a fresh object even when `orient` is unchanged, so an
    undiffed callback storms the parent during any resize.
  - `onRawSize` should likewise only fire when `w`/`h` actually change.
- `DetectorPanel`: forward both callbacks to `DetectorImage`, and pass the
  reported `orient` through to `DetectorRings` (the existing `orient` prop,
  currently never supplied). These are data callbacks / a render hint, not
  appearance — placement contract unaffected.

### 6. Focus page — `FocusPage.tsx`

- Hold detector `rawSize` (and `orient`) in state, set via the callbacks.
- Build `calibration = buildDetectorCalibration(experimentQ.data,
  rawSize)` and pass it to `DetectorPanel`.
- Derive and pass `imageAspect = rawSize ? rawSize.w / rawSize.h :
  undefined`. **Source it from the raw header dims, never from the decoded
  canvas/bitmap** — `resize_to_fit` floors to integers, so a bitmap-derived
  aspect drifts sub-pixel from the true raw aspect; keep this an explicit
  invariant. This is **required for registration**, not cosmetic: the SVG
  ring overlay stretches its `0 0 1 1` viewBox over the whole frame via
  `preserveAspectRatio="none"`. With the default square frame, a
  non-square image is letterboxed inside it (object-fit contain) while the
  rings stretch to the square — they would not register. Matching the
  frame aspect to the image makes the canvas fill the frame edge-to-edge
  so normalized ring coords map to image coords. (It also suppresses the
  landscape rotation, per Decisions.)
- **Two coupled first-paint transients, both benign and intended:**
  - Before headers arrive, `rawSize` is null → calibration null → centered
    fallback rings; once the header size arrives, the calibrated rings snap
    into place.
  - In the same pre-header window `imageAspect` is `undefined` → the frame
    is square → a non-square image is letterboxed and the *decorative*
    fallback rings are knowingly **not** registered to it (they stretch to
    the square frame). When `rawSize` lands the frame reshapes
    square→true-aspect (one reflow of the bounded frame box) and everything
    registers. The "registration is correct" guarantee applies only
    *after* headers arrive; the pre-header fallback rings are decorative by
    design, so this is acceptable.
  - There is also a ≤2-frame window where `evaluateOrient` can report
    `landscape` against the still-square frame (viewport ≥1400px, wide
    column) while calibration is still null; centered rings are
    rotation-invariant so it is invisible, and the `imageAspect` reshape
    then forces `portrait`. Confirm no flicker-rotate in the render-verify
    step at ≥1400px with a wide detector column.

### 7. Template + docs

- `packages/HimalayaUI/configs/simple.toml`: add commented example lines
  for the three new `[beamline]` fields.
- `docs/experiment-config.md`: document the three fields, the
  bottom-left pixel convention, and the "all four ingredients required"
  full-calibration rule in the `[beamline]` section.

## Data flow

```
experiment.toml [beamline]
  → config.jl (ExperimentConfig)              ← round-trip preserved
  → experiments.config blob (DB)
  → _beamline_from_config (routes_experiments) → experiment JSON
  → useExperiment() / api.ts Experiment ─┐
                                          ├→ buildDetectorCalibration(exp, rawSize)  ← pure, unit-tested
  image route X-Image-Width/Height ──────┘     (rawSize = raw detector px, from headers)
  (raw px, pre-resize)
  → DetectorPanel calibration + orient props
  → buildRingPlacements + qToImageRadius        (existing engine)
  → DetectorRings on the Focus page (rotates in lockstep via orient)
```

## Testing (TDD, six-layer contract)

- **config.jl** (`test/test_config.jl` or peer): round-trip
  `load_config` → `config_to_toml` → reparse preserves the three new
  fields; absent-fields case round-trips to `nothing`.
- **route** (`test/test_routes_experiments.jl` or peer):
  `_experiment_row_to_json` surfaces `beam_center_x/y` + `pixel_size_um`
  from a config blob; malformed/empty blob falls back without throwing.
- **image route** (`test/test_routes_image.jl`): the full-image response
  carries `X-Image-Width`/`X-Image-Height` equal to the raw raster size
  (assert against a known-size fixture larger than 1536px so the
  pre-resize capture is exercised, not just the no-op path).
- **frontend** (Vitest): `buildDetectorCalibration` — full calibration
  produced correctly (incl. µm→mm and m→mm conversions, and `imageSizePx`
  taken from the raw header size); each missing ingredient (and missing
  `rawSize`) → `null`; non-finite / zero / negative `rawSize` or numeric
  config → `null` (the finite/positive guard). Engine geometry is already
  covered by `detectorGeometry.test.ts`.
- **frontend regression**: the existing `DetectorImage.test.tsx` `fetch`
  mocks have **no `headers`** — they must grow a `headers: { get: () =>
  null }` stub (or the source's optional chaining keeps them green).
  Re-run the full `DetectorImage`/`DetectorPanel` suites; the new
  callbacks must not perturb the canvas-render / LUT / object-fit / orient
  cases.
- **live render verification**: serve a real experiment with a known beam
  center and confirm the rings center on the beamstop and a known-q ring
  lands on the corresponding diffraction radius. Confirms the bottom-left
  convention. **Two independent sources can produce a vertical-mirror
  symptom** — disambiguate them: (a) the bottom-left vs top-left beam
  convention (a one-line flip in `buildRingPlacements`), and (b) a
  beamline that writes bottom-up TIFFs so the *image itself* is flipped
  relative to the beam coordinate space. Also verify at viewport ≥1400px
  with a wide detector column that no flicker-rotate occurs on first paint.

## Out of scope

- Editing the beam center via the UI (config + reingest only, per the
  read-only experiment-directory contract).
- Per-exposure beam centers (beam center is a detector constant for an
  experiment).
- Beam-center calibration tooling / auto-detection from the image.
- Detector tilt / rotation (PONI rot1/rot2/rot3) — the existing engine
  assumes a normal-incidence flat detector.
```
