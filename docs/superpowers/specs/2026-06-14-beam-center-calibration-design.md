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
- **Image dimensions come from the loaded image**, not config or a new
  backend endpoint. The `size="full"` detector route returns the
  native-resolution PNG (no resize), so the displayed image's
  `naturalWidth/naturalHeight` equals the raw detector pixel space —
  the same coordinate space the beam center is expressed in.

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
  `flight_path_m`), so a round-trip preserves `nothing`.

### 2. Backend route — `packages/HimalayaUI/src/routes_experiments.jl`

- Generalize `_q_units_from_config(cfg_text) -> String` into a single
  blob extractor that returns `q_units`, `beam_center_x`,
  `beam_center_y`, and `pixel_size_um` from **one** `TOML.parse`.
- `_experiment_row_to_json` merges those four into the experiment JSON.
  `energy_kev` / `flight_path_m` already ride along via `SELECT *`.
- Keep the defensive contract: a non-string, empty, or malformed config
  must never 500 a list endpoint — missing/garbage values resolve to
  `q_units = "A-1"` and `nothing` for the numeric fields.

### 3. Frontend types — `packages/HimalayaUI/frontend/src/api.ts`

Add to `interface Experiment` (all `number | null`):
`beam_center_x`, `beam_center_y`, `pixel_size_um`, `energy_kev`,
`flight_path_m`. (The last two already flow in the JSON but are
currently untyped.)

### 4. Calibration assembly — `focusAdapters.ts` (new pure helper)

```ts
buildDetectorCalibration(
  experiment: Experiment | undefined,
  naturalSize: { w: number; h: number } | null,
): DetectorCalibration | null
```

- Returns `null` if `experiment`, `naturalSize`, or **any** of the four
  ingredients is missing/null.
- Unit conversions: `sampleDistanceMm = flight_path_m * 1000`,
  `pixelSizeMm = pixel_size_um / 1000`.
- Maps `beam_center_x/y` → `beamCenterPx`, `naturalSize` → `imageSizePx`,
  `energy_kev` → `energyKeV`.

This is the unit-tested seam; all arithmetic and null-guarding lives here,
not in the component.

### 5. Image dimensions — `DetectorImage` / `DetectorPanel`

- `DetectorImage`: add optional `onNaturalSize?: (w: number, h: number) =>
  void`, fired from the `<img>` `onLoad` using
  `naturalWidth/naturalHeight`.
- `DetectorPanel`: forward `onNaturalSize` to `DetectorImage` (placement
  contract unaffected; it is a data callback, not appearance).

### 6. Focus page — `FocusPage.tsx`

- Hold detector natural size in state, set via the `onNaturalSize`
  callback.
- Build `calibration = buildDetectorCalibration(experimentQ.data,
  naturalSize)` and pass it to `DetectorPanel`.
- Derive `imageAspect = naturalSize ? w / h : undefined` and pass it —
  a bonus fix for the current hardcoded-square assumption (non-square
  detectors render correctly).
- First paint (image not yet loaded) → `naturalSize` null → calibration
  null → centered fallback, then rings snap into place when the image
  loads. This is the intended, documented transient.

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
  → useExperiment() / api.ts Experiment
  → buildDetectorCalibration(exp, naturalSize) ← pure, unit-tested
  → DetectorPanel calibration prop
  → buildRingPlacements + qToImageRadius        (existing engine)
  → DetectorRings on the Focus page
```

## Testing (TDD, six-layer contract)

- **config.jl** (`test/test_config.jl` or peer): round-trip
  `load_config` → `config_to_toml` → reparse preserves the three new
  fields; absent-fields case round-trips to `nothing`.
- **route** (`test/test_routes_experiments.jl` or peer):
  `_experiment_row_to_json` surfaces `beam_center_x/y` + `pixel_size_um`
  from a config blob; malformed/empty blob falls back without throwing.
- **frontend** (Vitest): `buildDetectorCalibration` — full calibration
  produced correctly (incl. µm→mm and m→mm conversions); each missing
  ingredient (and missing `naturalSize`) → `null`. Engine geometry is
  already covered by `detectorGeometry.test.ts`.
- **live render verification**: serve a real experiment with a known beam
  center and confirm the rings center on the beamstop and a known-q ring
  lands on the corresponding diffraction radius. Confirms the bottom-left
  convention; reveals a vertical flip if the assumption is wrong.

## Out of scope

- Editing the beam center via the UI (config + reingest only, per the
  read-only experiment-directory contract).
- Per-exposure beam centers (beam center is a detector constant for an
  experiment).
- Beam-center calibration tooling / auto-detection from the image.
- Detector tilt / rotation (PONI rot1/rot2/rot3) — the existing engine
  assumes a normal-incidence flat detector.
```
