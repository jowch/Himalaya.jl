# Greenfield Detector Renderer (Batch 2a) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Commit trailer (every commit): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`

**Goal:** Build the greenfield detector renderer — a pure real-image canvas (`DetectorImage`) warmed by a Print-native colormap, and a *geometrically real* Debye–Scherrer ring overlay (`DetectorRings`) centered on the beam center and clipped to the frame — under a new design-guard-exempt `src/print/detector/` dir, honed in Storybook against real captured fixtures.

**Architecture:** Three pure pieces + two renderers. `detectorGeometry.ts` converts a SAXS calibration (beam center px, sample–detector distance, beam energy, pixel size) + a peak's q into a ring placement in **normalized image coordinates** (beam center 0–1, radius as a fraction of image width); with no calibration it falls back to the old presentational map. `detectorLut.ts` builds a 256-entry warm OKLCH colormap (black → umber → terracotta → gold → cream) — a reusable Print asset replacing the old two-endpoint linear ramp. `DetectorImage.tsx` paints the real image through that LUT. `DetectorRings.tsx` is a **dumb, resolution-independent** overlay: it takes a beam center + per-ring radii already in image coordinates, draws glow+sharp+hit circles, and lets the frame's `overflow:hidden` clip the arcs of an off-center beam. The physics lives only in the helper, unit-tested in one place; when real calibration data lands later, only the helper's inputs get wired.

**Why this shape:** The current overlay is uncalibrated **by design** — `src/lib/qRing.ts` states "There is NO detector geometry… the rings are PRESENTATIONAL… not a calibrated projection (#180 open decision 1)." It centers every ring at the viewBox middle `(50,50)`. Real rings are concentric about the *beam center*, which is generally off-center, so large-q rings exit the frame as partial arcs. This plan replaces the decorative overlay with one built around real geometry, parameterized now and gracefully degrading until the backend supplies calibration.

**Tech Stack:** React 18 + TypeScript strict (`exactOptionalPropertyTypes` ON) + Vite + Storybook CSF3 + Vitest/RTL under JSDOM. Plain-JSX SVG (the detector is not an x–y plot). Canvas 2D + `createImageBitmap` + `OffscreenCanvas` for the image; OKLCH→sRGB in code (no browser dependency, so the LUT is testable under JSDOM).

---

## Design-guard: a new exempt dir

`src/print/components/**` is enforced placement-only, but the detector renderer *must* emit appearance literals (colormap OKLCH stops, ring stroke colours, opacities, glow). That's the "rendering layer that paints pixels" category the guard already exempts (`figure-export/`, legacy detector/heatmap). The greenfield exempt set today is `components/ui/`, `print/ui/`, `print/plot/` (`scripts/check-design.mjs:136-142`). `print/plot/` is "traces only" and the detector isn't a trace, so we add a fourth prefix `print/detector/`.

## File structure

| File | Responsibility |
|---|---|
| `scripts/check-design.mjs` (modify) | Add `print/detector/` to `isExcluded()` |
| `test/check-design.test.ts` (modify) | Assert `print/detector/**` exempt + `print/components/**` still enforced |
| `src/print/detector/detectorLut.ts` (create) | 256-entry warm OKLCH→sRGB colormap (the Print detector LUT) |
| `src/print/detector/detectorGeometry.ts` (create) | Calibration + q → ring placements in normalized image coords; presentational fallback |
| `src/print/detector/DetectorImage.tsx` (create) | Canvas renderer: real image → Print LUT; orient rotation; thumb scroll-gate; placeholder |
| `src/print/detector/DetectorRings.tsx` (create) | Image-coord ring overlay: beam center + radii, glow+sharp+hit, clip off-frame, props hover |
| `src/print/detector/index.ts` (create) | Barrel |
| `src/print/detector/*.stories.tsx` (create) | LUT swatch, image (full/thumb/missing), ring vocabulary + off-center composed |
| `test/print-detector/*.test.tsx` / `*.test.ts` (create) | One test file per module |

**Out of scope (separate follow-on plans):** the tier-2 `DetectorPanel` composite that consumes this pair (builds the `/api/exposures/:id/image` URL, wires Zustand `hoveredQ`, derives the q-set from the assignment cart per `FocusDetectorPanel.tsx:52-77`, threads the real calibration into `detectorGeometry`, registers the overlay against the live (possibly rotated/letterboxed) canvas rect, renders the exposure switcher). And `CombChart`, `WaterfallChart`, `CardFigure`, `CleanFigure`, `CustomPreview`.

> **Registration scope note for this batch.** The renderer + helper are proven in Storybook with a hand-set beam center over the fixture image in a square box (the focus full-frame box is `aspect-square`). Pixel-exact registration against the *live* canvas when it rotates 90° / letterboxes is a composite concern (real calibration + the actual display transform), called out in the follow-on. This batch delivers: correct geometry math (helper), correct draw+clip (component), correct LUT (image).

## Conventions (Batch-1 / plot-engine recipe)

Per module (TDD): write the test → run → fail → implement → write the story (for components) → green on `npm test -- print-detector/<Name>` + `npm run lint:design` + `npx tsc --noEmit -p tsconfig.build.json` → commit. Renderers are **pure**: no Zustand, no Query, no URL-building inside. Tests assert via roles / `data-*` / text — never class strings. `exactOptionalPropertyTypes`: forward optional handlers via conditional spread; never pass explicit `undefined`.

---

## Task 1: Design-guard exemption for `src/print/detector/`

**Files:** Modify `scripts/check-design.mjs:136-142`; Test `test/check-design.test.ts`.

- [ ] **Step 1: Write the failing test** — append to `test/check-design.test.ts`:

```ts
describe("isExcluded via scanContent — print/detector authoring", () => {
  it("excludes src/print/detector from appearance rules", () => {
    expect(scanContent("print/detector/DetectorRings.tsx", `<circle stroke="oklch(0.9 0.045 80)" />`)).toEqual([]);
    expect(scanContent("print/detector/detectorLut.ts", `const stop = "oklch(0.55 0.12 45)";`)).toEqual([]);
  });
  it("still ENFORCES print/components (placement-only composites)", () => {
    expect(scanContent("print/components/Foo.tsx", `<span className="text-[11px]" />`).length).toBeGreaterThan(0);
  });
});
```

- [ ] **Step 2: Run → fail.** `npm test -- check-design` → FAIL (first case returns a violation; `print/detector/` not excluded yet).

- [ ] **Step 3: Add the prefix** in `scripts/check-design.mjs`:

```js
// src/components/ui/**, src/print/ui/**, src/print/plot/**, and src/print/detector/**
// are excluded (appearance authored — primitives, the trace-plot engine, and the
// detector rendering layer that paints real pixels).
function isExcluded(relPath) {
  return (
    relPath.startsWith("components/ui/") ||
    relPath.startsWith("print/ui/") ||
    relPath.startsWith("print/plot/") ||
    relPath.startsWith("print/detector/")
  );
}
```

- [ ] **Step 4: Run → pass.** `npm test -- check-design` → PASS (new + all pre-existing cases).

- [ ] **Step 5: Commit.** `git commit -m "Exempt src/print/detector/ from the design guard (rendering layer)"`

---

## Task 2: `detectorLut` — the Print detector colormap

**Files:** Create `src/print/detector/detectorLut.ts`; Test `test/print-detector/detectorLut.test.ts`.

A 256-entry RGB lookup table sampled from warm OKLCH stops, interpolated in OKLCH and converted to sRGB in code (Ottosson's transform). SAXS PNGs arrive log-normalized grayscale; the LUT maps intensity 0→warm near-black (window backing), rising through umber and a terracotta bloom to a cream highlight — a film-exposure feel, Print-native. Memoized (built once).

- [ ] **Step 1: Write the failing test** — `test/print-detector/detectorLut.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import { buildPrintDetectorLut, PRINT_DETECTOR_STOPS } from "../../src/print/detector/detectorLut";

describe("buildPrintDetectorLut", () => {
  const lut = buildPrintDetectorLut();

  it("has 256 RGB triples (768 bytes)", () => {
    expect(lut.length).toBe(256 * 3);
  });

  it("is non-inverting — luminance is monotonically non-decreasing", () => {
    const luma = (i: number) => 0.2126 * lut[i*3] + 0.7152 * lut[i*3+1] + 0.0722 * lut[i*3+2];
    for (let i = 1; i < 256; i++) {
      // Allow a tiny epsilon for rounding; the ramp must never step down.
      expect(luma(i)).toBeGreaterThanOrEqual(luma(i - 1) - 1.5);
    }
  });

  it("endpoints: index 0 is dark, index 255 is light and warm", () => {
    const sum0 = lut[0] + lut[1] + lut[2];
    const sum255 = lut[255*3] + lut[255*3+1] + lut[255*3+2];
    expect(sum0).toBeLessThan(180);          // near-black window backing
    expect(sum255).toBeGreaterThan(600);     // cream highlight
    expect(lut[255*3]).toBeGreaterThanOrEqual(lut[255*3+2]); // warm: R ≥ B
  });

  it("blooms warm in the mid-highs (terracotta: R clearly exceeds B near t≈0.62)", () => {
    const i = Math.round(0.62 * 255);
    expect(lut[i*3] - lut[i*3+2]).toBeGreaterThan(20);
  });

  it("exposes its stops for the Storybook swatch", () => {
    expect(PRINT_DETECTOR_STOPS.length).toBeGreaterThanOrEqual(4);
    expect(PRINT_DETECTOR_STOPS[0].t).toBe(0);
    expect(PRINT_DETECTOR_STOPS[PRINT_DETECTOR_STOPS.length - 1].t).toBe(1);
  });
});
```

- [ ] **Step 2: Run → fail.** `npm test -- print-detector/detectorLut` → FAIL (module missing).

- [ ] **Step 3: Implement `src/print/detector/detectorLut.ts`:**

```ts
// The Print detector colormap. SAXS images arrive log-normalized grayscale;
// this warm OKLCH ramp (black → umber → terracotta bloom → gold → cream) is the
// Print-native replacement for the old two-endpoint frame-edge→frame-signal lerp.
// OKLCH→sRGB is done in code (Ottosson) so the LUT is pure + JSDOM-testable.

export interface LutStop { t: number; L: number; C: number; h: number } // OKLCH

// Anchored to the Print warm axis (hues ~45–85° = terracotta→gold). Lightness
// rises monotonically so the ramp is non-inverting.
export const PRINT_DETECTOR_STOPS: readonly LutStop[] = [
  { t: 0.00, L: 0.150, C: 0.010, h: 55 }, // warm near-black (window backing)
  { t: 0.35, L: 0.320, C: 0.060, h: 50 }, // deep umber
  { t: 0.62, L: 0.550, C: 0.130, h: 45 }, // terracotta bloom (the Print accent heat)
  { t: 0.85, L: 0.780, C: 0.090, h: 70 }, // warm gold
  { t: 1.00, L: 0.930, C: 0.020, h: 85 }, // cream highlight
];

function clamp01(v: number): number { return v < 0 ? 0 : v > 1 ? 1 : v; }

// OKLCH → sRGB (0..255). Björn Ottosson's oklab matrices + sRGB gamma.
function oklchToSrgb(L: number, C: number, h: number): [number, number, number] {
  const hr = (h * Math.PI) / 180;
  const a = C * Math.cos(hr);
  const b = C * Math.sin(hr);
  const l_ = L + 0.3963377774 * a + 0.2158037573 * b;
  const m_ = L - 0.1055613458 * a - 0.0638541728 * b;
  const s_ = L - 0.0894841775 * a - 1.2914855480 * b;
  const l = l_ * l_ * l_, m = m_ * m_ * m_, s = s_ * s_ * s_;
  const lr =  4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s;
  const lg = -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s;
  const lb = -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s;
  const gamma = (v: number): number => {
    const x = v <= 0.0031308 ? 12.92 * v : 1.055 * Math.pow(v, 1 / 2.4) - 0.055;
    return clamp01(x);
  };
  return [
    Math.round(gamma(lr) * 255),
    Math.round(gamma(lg) * 255),
    Math.round(gamma(lb) * 255),
  ];
}

let cached: Uint8ClampedArray | null = null;

/** Build (and memoize) the 256×RGB Print detector LUT. */
export function buildPrintDetectorLut(): Uint8ClampedArray {
  if (cached) return cached;
  const lut = new Uint8ClampedArray(256 * 3);
  const stops = PRINT_DETECTOR_STOPS;
  for (let i = 0; i < 256; i++) {
    const t = i / 255;
    // Find the bracketing stop pair.
    let a = stops[0], b = stops[stops.length - 1];
    for (let k = 0; k < stops.length - 1; k++) {
      if (t >= stops[k].t && t <= stops[k + 1].t) { a = stops[k]; b = stops[k + 1]; break; }
    }
    const span = b.t - a.t || 1;
    const f = (t - a.t) / span;
    // Linear in L/C/h — the warm stops sit in a narrow hue band so no hue
    // wrap-around handling is needed.
    const [r, g, bl] = oklchToSrgb(
      a.L + f * (b.L - a.L),
      a.C + f * (b.C - a.C),
      a.h + f * (b.h - a.h),
    );
    lut[i * 3] = r; lut[i * 3 + 1] = g; lut[i * 3 + 2] = bl;
  }
  cached = lut;
  return lut;
}
```

- [ ] **Step 4: Run → pass** (`npm test -- print-detector/detectorLut`); `lint:design` clean; `tsc` clean.

- [ ] **Step 5: Commit.** `git commit -m "Print detector LUT: warm OKLCH colormap (black→terracotta→cream)"`

---

## Task 3: `detectorGeometry` — q → ring placements in image coordinates

**Files:** Create `src/print/detector/detectorGeometry.ts`; Test `test/print-detector/detectorGeometry.test.ts`.

Pure SAXS geometry. `qToImageRadius` does `λ = 12.398/E(keV)` → `θ = asin(qλ/4π)` → `r_mm = distance·tan(2θ)` → `r_px = r_mm/pixelSize` → normalize by image width. `buildRingPlacements` maps a q-set + optional calibration into `{ beamCenter, rings }` in normalized (0–1) image coordinates; with `null` calibration it falls back to the presentational map (centered, radius ∝ q-range) so the component works before real data exists.

- [ ] **Step 1: Write the failing test** — `test/print-detector/detectorGeometry.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import {
  qToImageRadius, buildRingPlacements, type DetectorCalibration,
} from "../../src/print/detector/detectorGeometry";

// λ = 12.398/12.398 = 1.0 Å; 1000 mm distance; 0.1 mm pixels; 1000×1000 image.
const CAL: DetectorCalibration = {
  beamCenterPx: { x: 500, y: 500 },
  imageSizePx: { w: 1000, h: 1000 },
  sampleDistanceMm: 1000,
  pixelSizeMm: 0.1,
  energyKeV: 12.398,
};

describe("qToImageRadius", () => {
  it("matches the hand-computed small-angle radius for q=0.1", () => {
    // θ = asin(0.1·1/(4π)) = 0.0079577 rad; 2θ → tan ≈ 0.015916;
    // r_mm = 1000·0.015916 = 15.916; r_px = 159.16; /1000 = 0.15916.
    expect(qToImageRadius(0.1, CAL)).toBeCloseTo(0.1592, 3);
  });
  it("grows monotonically with q", () => {
    expect(qToImageRadius(0.2, CAL)).toBeGreaterThan(qToImageRadius(0.1, CAL));
  });
});

describe("buildRingPlacements", () => {
  it("with calibration: beam center is the normalized px center; radii from geometry", () => {
    const { beamCenter, rings } = buildRingPlacements([0.1, 0.2], CAL);
    expect(beamCenter).toEqual({ x: 0.5, y: 0.5 });
    expect(rings.map((r) => r.q)).toEqual([0.1, 0.2]);
    expect(rings[0].r).toBeCloseTo(0.1592, 3);
  });
  it("off-center beam normalizes correctly", () => {
    const { beamCenter } = buildRingPlacements([0.1], { ...CAL, beamCenterPx: { x: 950, y: 60 } });
    expect(beamCenter.x).toBeCloseTo(0.95, 3);
    expect(beamCenter.y).toBeCloseTo(0.06, 3);
  });
  it("null calibration → presentational fallback: centered, radius ∝ q-range", () => {
    const { beamCenter, rings } = buildRingPlacements([0.1, 0.2, 0.3], null);
    expect(beamCenter).toEqual({ x: 0.5, y: 0.5 });
    // innermost ≈ RING_R_MIN/VIEWBOX, span ↑ with q; monotonic + within (0, 0.5].
    expect(rings[0].r).toBeGreaterThan(0);
    expect(rings[2].r).toBeGreaterThan(rings[0].r);
    expect(rings[2].r).toBeLessThanOrEqual(0.5);
  });
  it("carries color + ghost through from the q-set entries", () => {
    const { rings } = buildRingPlacements(
      [{ q: 0.1, color: "var(--color-success)" }, { q: 0.2, ghost: true }], CAL,
    );
    expect(rings[0].color).toBe("var(--color-success)");
    expect(rings[1].ghost).toBe(true);
  });
});
```

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement `src/print/detector/detectorGeometry.ts`:**

```ts
import { RING_R_MIN, RING_R_SPAN, RING_VIEWBOX } from "../../lib/qRing";

/** Standard SAXS/WAXS calibration (PONI-style). All lengths share one unit
 *  system; beamCenter + imageSize are in the displayed image's pixel space. */
export interface DetectorCalibration {
  beamCenterPx: { x: number; y: number };
  imageSizePx: { w: number; h: number };
  sampleDistanceMm: number;
  pixelSizeMm: number;
  energyKeV: number;
}

/** A ring to draw, in NORMALIZED image coordinates: radius is a fraction of the
 *  image width, so the component is resolution-independent. */
export interface RingPlacement {
  q: number;
  r: number;
  color?: string;
  ghost?: boolean;
}

/** A q-value plus optional render hints. A bare number is sugar for { q }. */
export type RingInput = number | { q: number; color?: string; ghost?: boolean };

const PLANCK_KEV_ANGSTROM = 12.39842; // hc in keV·Å

/** q (Å⁻¹) → ring radius as a fraction of image width, via SAXS geometry. */
export function qToImageRadius(q: number, cal: DetectorCalibration): number {
  const lambda = PLANCK_KEV_ANGSTROM / cal.energyKeV;          // Å
  const sinTheta = Math.min(1, (q * lambda) / (4 * Math.PI));
  const theta = Math.asin(sinTheta);
  const rMm = cal.sampleDistanceMm * Math.tan(2 * theta);
  const rPx = rMm / cal.pixelSizeMm;
  return rPx / cal.imageSizePx.w;
}

function norm(input: RingInput): { q: number; color?: string; ghost?: boolean } {
  return typeof input === "number" ? { q: input } : input;
}

/** Build ring placements (normalized image coords) + the beam center.
 *  `null` calibration → presentational fallback (centered, radius ∝ q-range),
 *  matching the legacy behaviour until real geometry data is wired. */
export function buildRingPlacements(
  inputs: RingInput[],
  cal: DetectorCalibration | null,
): { beamCenter: { x: number; y: number }; rings: RingPlacement[] } {
  const items = inputs.map(norm);
  if (cal) {
    return {
      beamCenter: {
        x: cal.beamCenterPx.x / cal.imageSizePx.w,
        y: cal.beamCenterPx.y / cal.imageSizePx.h,
      },
      rings: items.map((it) => ({
        q: it.q,
        r: qToImageRadius(it.q, cal),
        ...(it.color !== undefined ? { color: it.color } : {}),
        ...(it.ghost !== undefined ? { ghost: it.ghost } : {}),
      })),
    };
  }
  // Presentational fallback: normalize the old 100-unit viewBox radii to 0..1.
  const qs = items.map((it) => it.q);
  const qLo = qs.length ? Math.min(...qs) : 0;
  const qHi = qs.length ? Math.max(...qs) : 1;
  const span = qHi - qLo || 1;
  return {
    beamCenter: { x: 0.5, y: 0.5 },
    rings: items.map((it) => {
      const t = (it.q - qLo) / span;
      const rViewbox = RING_R_MIN + t * RING_R_SPAN; // 12..45 in the 100 box
      return {
        q: it.q,
        r: rViewbox / RING_VIEWBOX, // → 0.12..0.45 normalized
        ...(it.color !== undefined ? { color: it.color } : {}),
        ...(it.ghost !== undefined ? { ghost: it.ghost } : {}),
      };
    }),
  };
}
```

- [ ] **Step 4: Run → pass**; `lint:design` clean; `tsc` clean.

- [ ] **Step 5: Commit.** `git commit -m "detectorGeometry: SAXS q→image-coord ring placements + presentational fallback"`

---

## Task 4: `DetectorImage` — real-image canvas through the Print LUT

**Files:** Create `src/print/detector/DetectorImage.tsx`; Test `test/print-detector/DetectorImage.test.tsx`.

Ports `src/components/DetectorImage.tsx` with two changes: (1) the URL triple `{exposureId, imagePath, imageVersion}` collapses to one `src: string | null` (caller builds the URL → Storybook feeds a fixture asset; the future composite feeds the API URL); (2) the per-pixel warm ramp uses `buildPrintDetectorLut()` (index by intensity) instead of lerping two `getCssColor` endpoints. Orient rotation (`decideOrient`), IntersectionObserver thumb-gating, and the `frame-window` placeholder port unchanged.

- [ ] **Step 1: Write the failing tests** — `test/print-detector/DetectorImage.test.tsx`:

```tsx
import { render, screen, waitFor } from "@testing-library/react";
import { vi, beforeEach, test, expect } from "vitest";
import { DetectorImage } from "../../src/print/detector/DetectorImage";

const TINY_PNG =
  "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwADhQGAWjR9awAAAABJRU5ErkJggg==";

beforeEach(() => {
  global.fetch = vi.fn().mockResolvedValue({
    ok: true,
    blob: () => Promise.resolve(new Blob(
      [Uint8Array.from(atob(TINY_PNG), (c) => c.charCodeAt(0))], { type: "image/png" })),
  } as Response);
  global.createImageBitmap = vi.fn().mockResolvedValue({ width: 1, height: 1, close: vi.fn() } as unknown as ImageBitmap);
  const mockOffscreen = {
    getContext: () => ({ drawImage: vi.fn(), getImageData: () => ({ data: new Uint8ClampedArray(4) }) }),
  };
  // @ts-expect-error JSDOM stub
  global.OffscreenCanvas = vi.fn().mockImplementation(() => mockOffscreen);
});

test("shows the frame-window placeholder when src is null", () => {
  render(<DetectorImage src={null} size="full" />);
  expect(screen.getByTestId("detector-image-placeholder")).toHaveAttribute("data-variant", "frame-window");
});

test("renders a canvas (role=img) when src is provided", async () => {
  render(<DetectorImage src="/fixtures/thumbs/37.png" size="full" />);
  await waitFor(() => expect(screen.getByRole("img", { hidden: true })).toBeInTheDocument());
});

test("fetches exactly the src it was given (no URL building inside)", async () => {
  const spy = vi.fn().mockResolvedValue({
    ok: true, blob: () => Promise.resolve(new Blob([new Uint8Array(0)], { type: "image/png" })),
  } as Response);
  global.fetch = spy;
  render(<DetectorImage src="/api/exposures/42/image?thumb=1&v=v1-9" size="thumb" />);
  await waitFor(() => expect(spy).toHaveBeenCalled());
  expect(spy.mock.calls[0][0]).toBe("/api/exposures/42/image?thumb=1&v=v1-9");
  expect(spy.mock.calls[0][1]).toBeUndefined();
});

test("LUT is non-inverting — brighter source pixel → lighter output", async () => {
  // Two source pixels: intensity 0 and intensity 255 (R channel is read as t).
  const srcData = new Uint8ClampedArray([0,0,0,255, 255,255,255,255]);
  const mockOffscreen = { getContext: () => ({ drawImage: vi.fn(), getImageData: () => ({ data: srcData }) }) };
  // @ts-expect-error JSDOM stub
  global.OffscreenCanvas = vi.fn().mockImplementation(() => mockOffscreen);
  global.createImageBitmap = vi.fn().mockResolvedValue({ width: 2, height: 1, close: vi.fn() } as unknown as ImageBitmap);

  let captured: Uint8ClampedArray | null = null;
  vi.spyOn(HTMLCanvasElement.prototype, "getContext").mockImplementation(function () {
    return {
      drawImage: () => {},
      putImageData: (img: ImageData) => { captured = img.data; },
    } as unknown as CanvasRenderingContext2D;
  } as typeof HTMLCanvasElement.prototype.getContext);

  try {
    render(<DetectorImage src="/x.png" size="full" />);
    await waitFor(() => expect(captured).not.toBeNull());
    const d = captured as unknown as Uint8ClampedArray;
    const dark = d[0] + d[1] + d[2];
    const bright = d[4] + d[5] + d[6];
    expect(bright).toBeGreaterThan(dark);           // non-inverting
    expect(d[0]).toBeGreaterThanOrEqual(d[2]);       // warm at the dark end (R ≥ B)
  } finally {
    vi.restoreAllMocks();
  }
});
```

Then copy the U-3 thumb-locks-portrait pair and `forceWideGeometry` helper verbatim from `test/DetectorImage.test.tsx:119-199`, replacing each `<DetectorImage exposureId={N} imagePath="…" imageVersion="…" size="…" />` with `<DetectorImage src="/x.png" size="…" />`.

- [ ] **Step 2: Run → fail** (module missing).

- [ ] **Step 3: Implement `src/print/detector/DetectorImage.tsx`:**

```tsx
import { useCallback, useEffect, useRef, useState, type CSSProperties } from "react";
import { decideOrient } from "../../lib/detectorOrient";
import { buildPrintDetectorLut } from "./detectorLut";

interface Props {
  /** Image URL to fetch + paint; null → frame-window placeholder. Caller builds
   *  it (Storybook: a fixture asset URL; composite: the API image URL). */
  src: string | null;
  /** `thumb` locks portrait + scroll-gates the fetch; `full` rotates + eager. */
  size: "thumb" | "full";
  className?: string;
}

interface Layout { orient: "portrait" | "landscape"; caps: { maxW: number; maxH: number } | null }

export function DetectorImage({ src, size, className }: Props): JSX.Element {
  const wrapperRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [layout, setLayout] = useState<Layout>({ orient: "portrait", caps: null });
  const [hasIntersected, setHasIntersected] = useState<boolean>(
    () => size === "full" || typeof window === "undefined" ||
      typeof window.IntersectionObserver !== "function",
  );

  const evaluateOrient = useCallback(() => {
    const wrapper = wrapperRef.current, canvas = canvasRef.current;
    if (!wrapper || !canvas || !canvas.width || !canvas.height) return;
    if (size === "thumb") { setLayout({ orient: "portrait", caps: null }); return; }
    const cw = wrapper.clientWidth, ch = wrapper.clientHeight;
    if (cw === 0 || ch === 0) return;
    setLayout(decideOrient({
      containerW: cw, containerH: ch, imageW: canvas.width, imageH: canvas.height,
      viewportW: typeof window !== "undefined" ? window.innerWidth : 0,
    }));
  }, [size]);

  const renderImage = useCallback(async () => {
    const canvas = canvasRef.current;
    if (!canvas || !src || !hasIntersected) return;
    const res = await fetch(src);
    if (!res.ok) return;
    const bitmap = await createImageBitmap(await res.blob());
    const { width, height } = bitmap;
    const off = new OffscreenCanvas(width, height);
    const offCtx = off.getContext("2d")!;
    offCtx.drawImage(bitmap, 0, 0);
    bitmap.close();
    const imageData = offCtx.getImageData(0, 0, width, height);

    // Warm, NON-inverting Print LUT: index the 256-entry colormap by the
    // grayscale intensity (R channel of the log-normalized PNG).
    const lut = buildPrintDetectorLut();
    const data = imageData.data;
    for (let i = 0; i < data.length; i += 4) {
      const t = data[i];           // 0..255 intensity
      data[i]     = lut[t * 3];
      data[i + 1] = lut[t * 3 + 1];
      data[i + 2] = lut[t * 3 + 2];
      data[i + 3] = 255;
    }
    canvas.width = width;
    canvas.height = height;
    canvas.getContext("2d")?.putImageData(imageData, 0, 0);
    evaluateOrient();
  }, [src, hasIntersected, evaluateOrient]);

  useEffect(() => { renderImage(); }, [renderImage]);

  useEffect(() => {
    if (hasIntersected) return;
    const wrapper = wrapperRef.current;
    if (!wrapper) return;
    if (typeof window === "undefined" || typeof window.IntersectionObserver !== "function") {
      setHasIntersected(true); return;
    }
    const io = new window.IntersectionObserver((entries) => {
      for (const e of entries) if (e.isIntersecting) { setHasIntersected(true); io.disconnect(); break; }
    }, { rootMargin: "200px 0px" });
    io.observe(wrapper);
    return () => io.disconnect();
  }, [hasIntersected]);

  useEffect(() => {
    const wrapper = wrapperRef.current;
    if (!wrapper) return;
    const ro = new ResizeObserver(() => evaluateOrient());
    ro.observe(wrapper);
    return () => ro.disconnect();
  }, [evaluateOrient]);

  if (!src) {
    return (
      <div data-testid="detector-image-placeholder" data-variant="frame-window"
        className={`flex items-center justify-center bg-frame-edge text-frame-tag font-mono text-sm tracking-wide ${className ?? ""}`}>
        No image
      </div>
    );
  }

  const canvasStyle: CSSProperties = {
    imageRendering: "pixelated",
    ...(layout.orient === "landscape" && layout.caps
      ? { maxWidth: `${layout.caps.maxW}px`, maxHeight: `${layout.caps.maxH}px`,
          transform: "rotate(90deg)", transformOrigin: "center" }
      : { maxWidth: "100%", maxHeight: "100%" }),
  };

  return (
    <div ref={wrapperRef} data-orient={layout.orient}
         className={`flex items-center justify-center w-full h-full overflow-hidden ${className ?? ""}`}>
      <canvas ref={canvasRef} role="img" aria-label="Detector image" style={canvasStyle} />
    </div>
  );
}
```

> **Fixture coarseness (not a code change).** The fixture thumbs are 121×128 8-bit grayscale, so the LUT applies cleanly (no preprocessing needed to apply it — that resolves your "process them to apply a LUT" worry: grayscale-intensity is all the LUT needs). At a 420px `full` frame they upscale ~3.5×, so `imageRendering: pixelated` looks blocky — that's a fixture artifact only; production serves the real full-res image. We keep `pixelated` (preserves real pixel data — scientific honesty). If you want the Storybook full-frame to look smoother, that's a story-wrapper CSS choice, not a component change.

- [ ] **Step 4: Run → pass**; `lint:design` clean; `tsc` clean.

- [ ] **Step 5: Commit.** `git commit -m "Greenfield DetectorImage: pure canvas renderer through the Print LUT (src prop)"`

---

## Task 5: `DetectorImage` + LUT stories

**Files:** Create `src/print/detector/DetectorImage.stories.tsx`.

```tsx
import type { Meta, StoryObj } from "@storybook/react";
import { DetectorImage } from "./DetectorImage";
import { buildPrintDetectorLut } from "./detectorLut";
import thumb37 from "../fixtures/thumbs/37.png?url";

const meta: Meta<typeof DetectorImage> = { title: "detector/DetectorImage", component: DetectorImage };
export default meta;
type Story = StoryObj<typeof DetectorImage>;

const frame = "aspect-square w-[420px] overflow-hidden rounded border border-frame-edge bg-frame-edge";

export const Full: Story = {
  render: () => <div className={frame}><DetectorImage src={thumb37} size="full" className="h-full w-full" /></div>,
};
export const Thumb: Story = {
  render: () => (
    <div className="h-8 w-8 overflow-hidden rounded-sm border border-frame-edge bg-frame-edge">
      <DetectorImage src={thumb37} size="thumb" className="h-full w-full" />
    </div>
  ),
};
export const MissingImage: Story = {
  render: () => <div className={frame}><DetectorImage src={null} size="full" className="h-full w-full" /></div>,
};

// The colormap itself, as a horizontal swatch — for fidelity review of the ramp.
export const Colormap: Story = {
  render: () => {
    const lut = buildPrintDetectorLut();
    return (
      <div className="flex h-10 w-[512px] overflow-hidden rounded border border-hair">
        {Array.from({ length: 256 }, (_, i) => (
          <div key={i} className="h-full flex-1"
               style={{ background: `rgb(${lut[i*3]}, ${lut[i*3+1]}, ${lut[i*3+2]})` }} />
        ))}
      </div>
    );
  },
};
```

- [ ] **Verify:** `npm run build-storybook` → exit 0. (Manual: `npm run storybook` → `detector/DetectorImage/Colormap` shows the warm ramp; `Full` shows the real diffraction image in Print colour.) **Commit.**

---

## Task 6: `DetectorRings` — image-coordinate ring overlay

**Files:** Create `src/print/detector/DetectorRings.tsx`; Test `test/print-detector/DetectorRings.test.tsx`.

Dumb + resolution-independent: takes `beamCenter` (0–1) and `rings` (radii as fraction of image width, from `detectorGeometry`), draws each as glow+sharp+hit circles in a `viewBox="0 0 1 1"` SVG that overlays the image box, and lets the frame's `overflow:hidden` clip off-frame arcs. Hover is props-driven (`hoveredQ`/`onHoverQ`) — no Zustand. `hot` (q matches `hoveredQ`) → accent; `ghost` → dashed hollow; no colour → neutral ink-faint.

- [ ] **Step 1: Write the failing tests** — `test/print-detector/DetectorRings.test.tsx`:

```tsx
import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { DetectorRings } from "../../src/print/detector/DetectorRings";
import type { RingPlacement } from "../../src/print/detector/detectorGeometry";

const rings: RingPlacement[] = [
  { q: 0.10, r: 0.12, color: "var(--color-success)" },
  { q: 0.20, r: 0.22, color: "var(--color-success)" },
  { q: 0.30, r: 0.34, color: "var(--color-success)", ghost: true },
  { q: 0.40, r: 0.46 }, // leftover, neutral
];
const beamCenter = { x: 0.5, y: 0.5 };

describe("DetectorRings", () => {
  it("renders one ring group per placement, each with glow + sharp + hit", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={rings} />);
    expect(container.querySelectorAll('[data-role="det-ring"]').length).toBe(4);
    const g = container.querySelector('[data-role="det-ring"]')!;
    expect(g.querySelector('[data-role="ring-glow"]')).toBeTruthy();
    expect(g.querySelector('[data-role="ring-sharp"]')).toBeTruthy();
    expect(g.querySelector('[data-role="ring-hit"]')).toBeTruthy();
  });

  it("centers every ring at the beam center (not the viewBox middle)", () => {
    const { container } = render(
      <DetectorRings beamCenter={{ x: 0.8, y: 0.1 }} rings={[{ q: 0.1, r: 0.3 }]} />,
    );
    const sharp = container.querySelector('[data-role="ring-sharp"]')!;
    expect(Number(sharp.getAttribute("cx"))).toBeCloseTo(0.8, 5);
    expect(Number(sharp.getAttribute("cy"))).toBeCloseTo(0.1, 5);
    expect(Number(sharp.getAttribute("r"))).toBeCloseTo(0.3, 5);
  });

  it("draws a ghost ring dashed + hollow", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={rings} />);
    const ghost = container.querySelector('[data-ghost="true"] [data-role="ring-sharp"]')!;
    expect(ghost.getAttribute("stroke-dasharray")).toBeTruthy();
    expect(ghost.getAttribute("fill")).toBe("none");
  });

  it("a leftover ring (no color) strokes neutral ink-faint", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={[{ q: 0.4, r: 0.4 }]} />);
    expect(container.querySelector('[data-role="ring-sharp"]')!.getAttribute("stroke")).toBe("var(--color-ink-faint)");
  });

  it("the ring matching hoveredQ goes hot (accent stroke)", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={rings} hoveredQ={0.20} />);
    const hot = container.querySelector('[data-hot="true"]')!;
    expect(hot.getAttribute("data-ring-q")).toBe("0.2");
    expect(hot.querySelector('[data-role="ring-sharp"]')!.getAttribute("stroke")).toBe("var(--color-accent)");
  });

  it("no ring is hot when hoveredQ is undefined", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={rings} />);
    expect(container.querySelector('[data-hot="true"]')).toBeNull();
  });

  it("calls onHoverQ(q) on hit enter and onHoverQ(undefined) on leave", () => {
    const spy = vi.fn();
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={rings} onHoverQ={spy} />);
    const hit = container.querySelectorAll('[data-role="ring-hit"]')[1]; // q=0.20
    fireEvent.mouseEnter(hit); expect(spy).toHaveBeenCalledWith(0.20);
    fireEvent.mouseLeave(hit); expect(spy).toHaveBeenCalledWith(undefined);
  });

  it("hit ring is inert (pointer-events none) when onHoverQ is absent", () => {
    const { container } = render(<DetectorRings beamCenter={beamCenter} rings={rings} />);
    const hit = container.querySelector('[data-role="ring-hit"]') as SVGElement;
    expect(hit.style.pointerEvents).toBe("none");
  });
});
```

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement `src/print/detector/DetectorRings.tsx`:**

```tsx
import { type CSSProperties } from "react";
import type { RingPlacement } from "./detectorGeometry";

interface Props {
  /** Beam center in normalized image coords (0–1). The ring origin. */
  beamCenter: { x: number; y: number };
  /** Rings with radii in normalized image coords (fraction of image width). */
  rings: RingPlacement[];
  /** q hovered elsewhere (trace peak / comb). The matching ring lights to the
   *  accent q-link. Props-driven so the renderer is pure; the composite threads
   *  Zustand `hoveredQ`. */
  hoveredQ?: number;
  /** Fired on hit-ring enter (q) / leave (undefined). Omit → inert rings. */
  onHoverQ?: (q?: number) => void;
  /** Orientation when the parent drives it (mirrors DetectorImage). */
  orient?: "portrait" | "landscape";
}

export function DetectorRings({ beamCenter, rings, hoveredQ, onHoverQ, orient }: Props): JSX.Element {
  const svgStyle: CSSProperties = {
    position: "absolute", inset: 0, width: "100%", height: "100%",
    pointerEvents: "none", // hit rings re-enable individually
    ...(orient === "landscape" ? { transform: "rotate(90deg)", transformOrigin: "center" } : {}),
  };
  const TOL = 1e-6;

  return (
    <div className="pointer-events-none absolute inset-0 z-10">
      {/* viewBox is the normalized image box; preserveAspectRatio="none" stretches
          it to the displayed image rect so beam center + radii register. The frame's
          overflow:hidden clips the arcs of an off-center beam. */}
      <svg data-testid="detector-rings" data-orient={orient ?? "portrait"}
           viewBox="0 0 1 1" preserveAspectRatio="none" style={svgStyle} aria-hidden="true">
        {rings.map(({ q, r, color, ghost }) => {
          const hot = hoveredQ !== undefined && Math.abs(q - hoveredQ) <= TOL;
          const stroke = hot ? "var(--color-accent)" : (color ?? "var(--color-ink-faint)");
          // Radii/strokes are in the normalized (0..1) space; vector-effect keeps
          // stroke widths crisp despite preserveAspectRatio="none" scaling.
          const sharpW = hot ? 1.8 : ghost ? 0.8 : 1.0;
          const sharpOp = hot ? 0.95 : ghost ? 0.45 : 0.7;
          return (
            <g key={q} data-role="det-ring" data-ring-q={q}
               data-hot={hot ? "true" : undefined} data-ghost={ghost ? "true" : undefined}>
              <circle data-role="ring-glow" cx={beamCenter.x} cy={beamCenter.y} r={r}
                fill="none" stroke={stroke} strokeWidth={2.4} vectorEffect="non-scaling-stroke"
                opacity={hot ? 0.4 : ghost ? 0.06 : 0.18} style={{ pointerEvents: "none" }} />
              <circle data-role="ring-sharp" cx={beamCenter.x} cy={beamCenter.y} r={r}
                fill="none" stroke={stroke} strokeWidth={sharpW} vectorEffect="non-scaling-stroke"
                strokeDasharray={ghost ? "2 2.5" : undefined} opacity={sharpOp}
                style={{ pointerEvents: "none" }} />
              <circle data-role="ring-hit" cx={beamCenter.x} cy={beamCenter.y} r={r}
                fill="none" stroke="transparent" strokeWidth={5} vectorEffect="non-scaling-stroke"
                style={{ pointerEvents: onHoverQ ? "stroke" : "none", cursor: onHoverQ ? "pointer" : "default" }}
                {...(onHoverQ ? { onMouseEnter: () => onHoverQ(q), onMouseLeave: () => onHoverQ(undefined) } : {})} />
            </g>
          );
        })}
      </svg>
    </div>
  );
}
```

> The hit-ring's leave emits a raw `onHoverQ(undefined)`; the composite's `setHoveredQ` wrapper absorbs any A→B crossing transient (the legacy store-read guard does not belong in a pure renderer).

- [ ] **Step 4: Run → pass**; `lint:design` clean; `tsc` clean.

- [ ] **Step 5: Commit.** `git commit -m "Greenfield DetectorRings: image-coord overlay (beam center, off-frame clip, accent q-link)"`

---

## Task 7: `DetectorRings` story + composed off-center detector

**Files:** Create `src/print/detector/DetectorRings.stories.tsx`.

```tsx
import type { Meta, StoryObj } from "@storybook/react";
import { useState } from "react";
import { DetectorRings } from "./DetectorRings";
import { DetectorImage } from "./DetectorImage";
import { buildRingPlacements } from "./detectorGeometry";
import thumb37 from "../fixtures/thumbs/37.png?url";

const meta: Meta<typeof DetectorRings> = { title: "detector/DetectorRings", component: DetectorRings };
export default meta;
type Story = StoryObj<typeof DetectorRings>;

const frame = "relative aspect-square w-[420px] overflow-hidden rounded border border-frame-edge bg-frame-edge";
const QS = [
  { q: 0.082, color: "var(--color-success)" }, { q: 0.116, color: "var(--color-success)" },
  { q: 0.142, color: "var(--color-success)" }, { q: 0.171, color: "var(--color-success)", ghost: true },
  { q: 0.205 },
];

// Centered presentational fallback (null calibration) over the dark window.
export const Vocabulary: Story = {
  render: () => {
    const { beamCenter, rings } = buildRingPlacements(QS, null);
    return <div className={frame}><DetectorRings beamCenter={beamCenter} rings={rings} /></div>;
  },
};

// Off-center beam → partial arcs clipped by the frame. Demonstrates the geometry
// path with a hand-set calibration (beam near the bottom-right).
export const OffCenterBeam: Story = {
  render: () => {
    const { beamCenter, rings } = buildRingPlacements(QS, {
      beamCenterPx: { x: 880, y: 980 }, imageSizePx: { w: 1000, h: 1000 },
      sampleDistanceMm: 800, pixelSizeMm: 0.172, energyKeV: 12.4,
    });
    return <div className={frame}><DetectorRings beamCenter={beamCenter} rings={rings} /></div>;
  },
};

// The real composition: image + rings + interactive hover.
export const OverImage: Story = {
  render: () => {
    const [hovered, setHovered] = useState<number | undefined>(undefined);
    const { beamCenter, rings } = buildRingPlacements(QS, null);
    return (
      <div className={frame}>
        <DetectorImage src={thumb37} size="full" className="h-full w-full" />
        <DetectorRings beamCenter={beamCenter} rings={rings} hoveredQ={hovered} onHoverQ={setHovered} />
      </div>
    );
  },
};
```

- [ ] **Verify:** `npm run build-storybook` → exit 0. (Manual: `OffCenterBeam` shows partial arcs clipped at the frame edges; `OverImage` lights a ring terracotta on hover.) **Commit.**

---

## Task 8: Barrel + batch gates

**Files:** Create `src/print/detector/index.ts`.

```ts
export { DetectorImage } from "./DetectorImage";
export { DetectorRings } from "./DetectorRings";
export { buildPrintDetectorLut, PRINT_DETECTOR_STOPS } from "./detectorLut";
export {
  buildRingPlacements, qToImageRadius,
  type DetectorCalibration, type RingPlacement, type RingInput,
} from "./detectorGeometry";
```

- [ ] **Run the batch gates:**

```bash
cd /Users/me/projects/Himalaya.jl/.claude/worktrees/greenfield-ui-rebuild/packages/HimalayaUI/frontend
npm test -- print-detector       # all four test files green
npm run lint:design              # clean (detector/ exempt, components/ still enforced)
npm run build                    # tsc --noEmit + vite build, exit 0
npm run build-storybook          # exit 0
```

- [ ] **Commit.** `git commit -m "Greenfield detector renderer: barrel + batch gates green"`

---

## Verification

Per task: `npm test -- print-detector/<Name>` green; `npm run lint:design` clean; `npx tsc --noEmit -p tsconfig.build.json` clean. After the batch: `npm run build` + `npm run build-storybook` exit 0. Manual fidelity in Storybook: the `Colormap` ramp; the real image in Print colour; `OffCenterBeam` partial arcs; `OverImage` hover q-link — against `docs/redesign-mockups/2026-05-29-focus-plot.html` (`renderDetector`, lines 811-845).

## Follow-on (next plans)

- **Tier-2 `DetectorPanel` composite** (`src/print/components/`): API image URL, Zustand `hoveredQ` wiring, q-set from the assignment cart, **real `DetectorCalibration` threaded into `buildRingPlacements`**, registration against the live (rotating/letterboxing) canvas rect, exposure switcher + "Set rep".
- `CombChart`, `WaterfallChart`, `CardFigure`, `CleanFigure`, `CustomPreview` — each its own renderer plan.
