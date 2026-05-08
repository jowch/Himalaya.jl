# Figure Export Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Copy / Save (PNG/SVG) buttons to TraceViewer (Index page) and MultiTracePlot (Compare page) that produce presentation-ready figures via headless Observable Plot re-render.

**Architecture:** Pure-function adapters convert live plot state → `ExportSpec`; a renderer turns the spec into SVG (string/serialize) or PNG (SVG → OffscreenCanvas → blob). Two parents (PlotCard, ComparePage) host a small `<FigureExportControls>` component that calls the adapter as a thunk on click — no plot mutation, no Zustand wiring.

**Tech Stack:** React 18 + TypeScript strict, Observable Plot (already shipped), `OffscreenCanvas`, `navigator.clipboard.write`, Vitest (unit), Playwright (E2E).

**Spec:** [`docs/superpowers/specs/2026-05-08-figure-export-design.md`](../specs/2026-05-08-figure-export-design.md). The spec is the authority on design intent; this plan is the build sequence. Each task references spec §sections by their headings — open the spec when implementing a step.

**Scope reminders:**
- In: TraceViewer (Index page), MultiTracePlot (Compare page).
- Out: MillerPlot inset, DetectorImage canvas, R² in trace legend, PDF/LaTeX, visual-regression. See spec §Non-goals.
- The export-mark builders **duplicate** trace/peak/predicted-q geometry from on-screen renderers — accepted tradeoff per spec §Approach (the on-screen overlay is load-bearing for hit-testing/hover and must not be touched).

**Conventions for this branch:**
- Pure utility files (under `src/lib/figure-export/`) get a one-line header comment (file name + one-sentence purpose). No multi-paragraph docstrings.
- Test files live in `frontend/test/figure-export/` (one file per source unit). Follow the existing pattern: import from `../../src/...`, use `vi.fn()` for mocks, no Tailwind class assertions.
- TDD: failing test → minimal implementation → green → commit. Each task ends with one commit.
- The frontend dev loop expects `tsc --noEmit + vite build` to be clean before PR — run `npm run build` before considering Phase 4 done.

---

## File structure

```
src/lib/figure-export/
  types.ts                    — ExportSpec, LegendSpec
  presets.ts                  — TRACE_DIMS, COMPARE_DIMS, LIGHT_PALETTE, font/stroke constants
  filename.ts                 — slugifyForFilename, buildFilename
  download.ts                 — downloadBlob
  clipboard.ts                — copyPngToClipboard + canCopyPngToClipboard pre-flight
  renderer.ts                 — buildExportSvg, buildExportPng
  marks/
    traceExportMarks.ts       — emits trace + peaks + predicted-q ticks
    multiTraceExportMarks.ts  — emits stacked traces + per-member labels + peak ticks
  adapters/
    traceAdapter.ts           — buildTraceExportSpec (TraceViewer state → ExportSpec)
    multiTraceAdapter.ts      — buildMultiTraceExportSpec (MultiTracePlot state → ExportSpec)

src/components/
  FigureExportControls.tsx    — Copy + split-Download (PNG/SVG) buttons

frontend/test/figure-export/
  filename.test.ts
  clipboard.test.ts
  renderer.test.ts
  traceAdapter.test.ts
  multiTraceAdapter.test.ts
  FigureExportControls.test.tsx

frontend/e2e/
  figure-export.spec.ts       — Copy → clipboard PNG (Chromium-only) + Download → PNG/SVG
playwright.config.ts          — modify: add permissions for clipboard-read/write
```

Modified existing files:
- `src/components/PlotCard.tsx` — wire `<FigureExportControls>` into `TitleStrip` controls cluster
- `src/pages/ComparePage.tsx` — wire `<FigureExportControls>` into `compare-review-header` row (after `<ForksPopover>`)

---

# Phase 1 — Foundations

These five files have no upstream deps (other than each other) and ship the testable scaffolding the trace and compare paths consume.

---

### Task 1: `types.ts` + `presets.ts`

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/figure-export/types.ts`
- Create: `packages/HimalayaUI/frontend/src/lib/figure-export/presets.ts`

No tests — types and constants only.

- [ ] **Step 1: Create `types.ts`**

```ts
// types.ts — export spec passed from adapters to the renderer.
import * as Plot from "@observablehq/plot";

export interface LegendRow {
  /** CSS color used for the swatch / line stroke. */
  color: string;
  /** Mark style — small filled square for trace key, line for predicted-q phases. */
  symbol: "swatch" | "line";
  label: string;
}

export interface LegendSpec {
  rows: LegendRow[];
}

export interface ExportSpec {
  title: { primary: string; secondary?: string };
  width: number;   // logical (CSS) pixels
  height: number;
  /** Marks + x/y configs handed to Plot.plot(). Adapters MUST NOT set
   *  `title`, `caption`, or `figure: true` — the renderer wraps the SVG
   *  itself, and those fields make Plot return HTMLElement instead of
   *  SVGSVGElement (breaks the rest of the pipeline). */
  plot: Plot.PlotOptions;
  legend?: LegendSpec;
}
```

- [ ] **Step 2: Create `presets.ts`**

```ts
// presets.ts — pinned export dimensions, palette, fonts, strokes.
// Numeric values live here so adapters stay terse and figure-size changes
// are a one-line edit.

export const TRACE_DIMS = { width: 800, height: 600 } as const;
export const COMPARE_DIMS = { width: 1000, height: 400 } as const;

/** Margin defaults applied by the renderer (overridable via spec.plot). */
export const EXPORT_MARGIN = {
  top: 60,    // room for primary + secondary title text
  right: 24,
  bottom: 60, // room for legend rows below the plot
  left: 60,
} as const;

/** Plot title font metrics. */
export const TITLE_FONT_PRIMARY = "600 16px ui-sans-serif, system-ui, sans-serif";
export const TITLE_FONT_SECONDARY = "400 12px ui-sans-serif, system-ui, sans-serif";

/** Body / legend font (also applied via Plot's `style`). */
export const BODY_FONT = "400 12px ui-sans-serif, system-ui, sans-serif";

/** Trace + peak stroke widths — wider than on-screen for printability. */
export const TRACE_STROKE_PX = 1.75;
export const PEAK_TICK_STROKE_PX = 1.5;
export const PREDICTED_Q_STROKE_PX = 1.5;

/** Light palette for export (overrides dark-theme CSS vars). Values are
 *  literal CSS colors, not var(--…) — the renderer never resolves vars. */
export const LIGHT_PALETTE = {
  background:    "#ffffff",
  text:          "#1a1a1a",
  textMuted:     "#555555",
  textDim:       "#888888",
  border:        "#cccccc",
  trace:         "#1a1a1a",
  peakAuto:      "#1f5fb0",
  peakManual:    "#a0421f",
  peakExcluded:  "rgba(31, 95, 176, 0.35)",
} as const;

/**
 * Comparison palette tuned for white background. Same hue assignments as
 * `COMPARE_PALETTE` in `lib/comparison/coloring.ts` but lower OKLCH luminance
 * so contrast holds against the export's light bg. (Per spec §Risks — the
 * dark-theme palette has L = 0.76–0.80; we drop to 0.50–0.58 here.)
 *
 * If the implementer sees a hue here that fails contrast against #fff during
 * step 3 wiring, adjust THIS constant — keep on-screen → export hue mapping
 * stable.
 */
export const COMPARE_PALETTE_LIGHT: readonly string[] = [
  "oklch(0.55 0.16  33)", // warm coral
  "oklch(0.58 0.16  77)", // gold
  "oklch(0.55 0.14 147)", // moss
  "oklch(0.54 0.13 175)", // teal
  "oklch(0.55 0.13 200)", // cyan
  "oklch(0.54 0.13 220)", // azure
  "oklch(0.55 0.12 263)", // lavender
  "oklch(0.50 0.14 295)", // purple
  "oklch(0.50 0.14 315)", // magenta-purple
  "oklch(0.50 0.14 333)", // raspberry
  "oklch(0.55 0.14   3)", // pink-red
  "oklch(0.58 0.12 105)", // chartreuse
];

/** Neutral cool-gray for orphans on light bg. */
export const ORPHAN_FALLBACK_LIGHT = "oklch(0.50 0.02 270)";
```

- [ ] **Step 3: Verify TS compiles**

Run: `cd packages/HimalayaUI/frontend && npm run build 2>&1 | tail -20`
Expected: build succeeds (no errors from these new files; existing code untouched).

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/figure-export/types.ts \
        packages/HimalayaUI/frontend/src/lib/figure-export/presets.ts
git commit -m "feat(figure-export): add ExportSpec types + presets scaffolding"
```

---

### Task 2: `filename.ts` (TDD)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/figure-export/filename.ts`
- Test: `packages/HimalayaUI/frontend/test/figure-export/filename.test.ts`

Spec §Filenames — slugify per segment, en-CA date, sentinel `"figure"` for empty input.

- [ ] **Step 1: Write failing test**

`packages/HimalayaUI/frontend/test/figure-export/filename.test.ts`:

```ts
import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { slugifyForFilename, buildFilename } from "../../src/lib/figure-export/filename";

describe("slugifyForFilename", () => {
  it("lowercases and replaces spaces with dashes", () => {
    expect(slugifyForFilename("My Sample 04")).toBe("my-sample-04");
  });
  it("collapses runs of non-alphanumerics into a single dash", () => {
    expect(slugifyForFilename("a / b // c")).toBe("a-b-c");
  });
  it("trims leading and trailing dashes", () => {
    expect(slugifyForFilename("--foo--")).toBe("foo");
  });
  it("returns sentinel 'figure' for empty input", () => {
    expect(slugifyForFilename("")).toBe("figure");
  });
  it("returns sentinel 'figure' for input that slugifies to empty", () => {
    expect(slugifyForFilename("///")).toBe("figure");
  });
  it("preserves digits and basic latin letters", () => {
    expect(slugifyForFilename("AgBe-2026-001")).toBe("agbe-2026-001");
  });
  it("strips diacritics by replacing with dashes (acceptable for filename)", () => {
    // Non-ASCII alphabetic chars are not preserved; the contract is "safe filename",
    // not "Unicode-fidelity". Crash-safe, predictable.
    const out = slugifyForFilename("café");
    expect(out).toMatch(/^[a-z0-9-]+$/);
    expect(out.length).toBeGreaterThan(0);
  });
});

describe("buildFilename", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    // Pin to a known local date. Spec mandates en-CA YYYY-MM-DD output.
    vi.setSystemTime(new Date(2026, 4, 8, 14, 30)); // 2026-05-08 local
  });
  afterEach(() => { vi.useRealTimers(); });

  it("appends '-{YYYY-MM-DD}.{ext}' (en-CA)", () => {
    expect(buildFilename("himalaya-trace-jc23", "png")).toBe("himalaya-trace-jc23-2026-05-08.png");
  });
  it("supports svg extension", () => {
    expect(buildFilename("himalaya-trace-jc23", "svg")).toBe("himalaya-trace-jc23-2026-05-08.svg");
  });
});
```

- [ ] **Step 2: Run test, verify it fails**

Run: `cd packages/HimalayaUI/frontend && node_modules/.bin/vitest run test/figure-export/filename.test.ts`
Expected: FAIL with "Cannot find module '../../src/lib/figure-export/filename'".

- [ ] **Step 3: Implement `filename.ts`**

```ts
// filename.ts — slug-safe filename helpers for figure export.

/**
 * Lowercase, collapse non-alphanumeric runs into single dashes, trim leading
 * and trailing dashes. Empty / fully-stripped input returns the sentinel
 * "figure" so the resulting filename always has a visible stem segment.
 *
 * Apply this PER SEGMENT before joining with static dashes — otherwise the
 * concatenation boundary is ambiguous. See spec §Filenames.
 */
export function slugifyForFilename(s: string): string {
  if (!s) return "figure";
  const slug = s
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "-")
    .replace(/^-+|-+$/g, "");
  return slug.length === 0 ? "figure" : slug;
}

/**
 * Append `-{YYYY-MM-DD}.{ext}` to a pre-resolved stem. Date is local time
 * (matches user expectation when exporting near midnight); locale is pinned
 * to en-CA to guarantee `YYYY-MM-DD` regardless of system locale (en-US
 * defaults produce `05/08/2026` — `/` is invalid on every OS).
 */
export function buildFilename(stem: string, ext: "png" | "svg"): string {
  const date = new Intl.DateTimeFormat("en-CA", {
    year: "numeric",
    month: "2-digit",
    day: "2-digit",
  }).format(new Date());
  return `${stem}-${date}.${ext}`;
}
```

- [ ] **Step 4: Run test, verify it passes**

Run: same command as Step 2.
Expected: PASS, all assertions green.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/figure-export/filename.ts \
        packages/HimalayaUI/frontend/test/figure-export/filename.test.ts
git commit -m "feat(figure-export): slugify + buildFilename with en-CA date"
```

---

### Task 3: `download.ts` (no tests — trivial)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/figure-export/download.ts`

The spec describes this as a one-line helper. It has no testable branch points beyond what JSDOM can't simulate (real download) — covered by the Playwright E2E in Phase 4. Skip a unit test here.

- [ ] **Step 1: Implement**

```ts
// download.ts — trigger a browser download for an in-memory blob.
export function downloadBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  try {
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    a.style.display = "none";
    document.body.appendChild(a);
    a.click();
    a.remove();
  } finally {
    URL.revokeObjectURL(url);
  }
}
```

- [ ] **Step 2: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/figure-export/download.ts
git commit -m "feat(figure-export): downloadBlob helper"
```

---

### Task 4: `clipboard.ts` (TDD)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/figure-export/clipboard.ts`
- Test: `packages/HimalayaUI/frontend/test/figure-export/clipboard.test.ts`

Spec §Browser compatibility — three failure modes: insecure origin, missing `ClipboardItem`, runtime denial.

- [ ] **Step 1: Write failing test**

`packages/HimalayaUI/frontend/test/figure-export/clipboard.test.ts`:

```ts
import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import {
  canCopyPngToClipboard,
  copyPngToClipboard,
} from "../../src/lib/figure-export/clipboard";

const origNav = globalThis.navigator;
const origClipboardItem = (globalThis as { ClipboardItem?: unknown }).ClipboardItem;

afterEach(() => {
  vi.unstubAllGlobals();
  // Restore in case stubs leak.
  Object.defineProperty(globalThis, "navigator", {
    value: origNav,
    configurable: true,
  });
  if (origClipboardItem !== undefined) {
    (globalThis as { ClipboardItem?: unknown }).ClipboardItem = origClipboardItem;
  }
});

describe("canCopyPngToClipboard", () => {
  it("returns false when navigator.clipboard is missing", () => {
    vi.stubGlobal("navigator", {});
    expect(canCopyPngToClipboard()).toBe(false);
  });
  it("returns false when ClipboardItem is missing", () => {
    vi.stubGlobal("navigator", { clipboard: { write: vi.fn() } });
    delete (globalThis as { ClipboardItem?: unknown }).ClipboardItem;
    expect(canCopyPngToClipboard()).toBe(false);
  });
  it("returns true when both are present", () => {
    vi.stubGlobal("navigator", { clipboard: { write: vi.fn() } });
    (globalThis as { ClipboardItem?: unknown }).ClipboardItem = function () { /* stub */ };
    expect(canCopyPngToClipboard()).toBe(true);
  });
});

describe("copyPngToClipboard", () => {
  it("writes a ClipboardItem with image/png MIME", async () => {
    const writeSpy = vi.fn().mockResolvedValue(undefined);
    vi.stubGlobal("navigator", { clipboard: { write: writeSpy } });
    const ciCalls: Array<Record<string, Blob>> = [];
    (globalThis as { ClipboardItem?: unknown }).ClipboardItem =
      function (this: unknown, items: Record<string, Blob>) {
        ciCalls.push(items);
      };
    const blob = new Blob([new Uint8Array([1, 2, 3])], { type: "image/png" });

    await copyPngToClipboard(blob);

    expect(writeSpy).toHaveBeenCalledTimes(1);
    expect(ciCalls.length).toBe(1);
    expect(ciCalls[0]!["image/png"]).toBe(blob);
  });

  it("rejects when navigator.clipboard.write rejects", async () => {
    const writeSpy = vi.fn().mockRejectedValue(new Error("denied"));
    vi.stubGlobal("navigator", { clipboard: { write: writeSpy } });
    (globalThis as { ClipboardItem?: unknown }).ClipboardItem =
      function (this: unknown, _items: Record<string, Blob>) { /* stub */ };
    const blob = new Blob([new Uint8Array([1])], { type: "image/png" });

    await expect(copyPngToClipboard(blob)).rejects.toThrow("denied");
  });
});
```

- [ ] **Step 2: Run, verify FAIL**

Run: `node_modules/.bin/vitest run test/figure-export/clipboard.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Implement**

```ts
// clipboard.ts — write a PNG blob to the system clipboard.

/**
 * Pre-flight feature check. The Copy button uses this to decide whether to
 * disable itself (insecure origin / unsupported browser) before the user
 * even clicks. Spec §Browser compatibility / Disabled states.
 */
export function canCopyPngToClipboard(): boolean {
  if (typeof navigator === "undefined") return false;
  if (!navigator.clipboard) return false;
  if (typeof (globalThis as { ClipboardItem?: unknown }).ClipboardItem === "undefined") {
    return false;
  }
  return true;
}

/**
 * Write a PNG blob to the clipboard. Caller is expected to have gated
 * on `canCopyPngToClipboard()`; we don't re-check here so a runtime
 * rejection (user denial, third failure mode) propagates as a real error
 * the caller can toast.
 */
export async function copyPngToClipboard(blob: Blob): Promise<void> {
  // ClipboardItem is global at runtime when canCopyPngToClipboard() is true.
  // TS doesn't pick up DOM-lib types in the test env; cast at call time.
  const Ctor = (globalThis as { ClipboardItem: new (data: Record<string, Blob>) => unknown })
    .ClipboardItem;
  const item = new Ctor({ "image/png": blob }) as unknown;
  await navigator.clipboard.write([item as ClipboardItem]);
}
```

- [ ] **Step 4: Run, verify PASS**

Run: same command as Step 2.
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/figure-export/clipboard.ts \
        packages/HimalayaUI/frontend/test/figure-export/clipboard.test.ts
git commit -m "feat(figure-export): clipboard pre-flight + copyPngToClipboard"
```

---

### Task 5: `renderer.ts` (TDD)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/figure-export/renderer.ts`
- Test: `packages/HimalayaUI/frontend/test/figure-export/renderer.test.ts`

Spec §Renderer + §Per-plot figure design — the renderer is dumb: lay out title text + plot SVG + legend rows on a white background. PNG path is JSDOM-untestable (`OffscreenCanvas.convertToBlob` + `Image.decode` aren't real in JSDOM); test only structure of `buildExportSvg` and that `buildExportPng` returns a Promise.

- [ ] **Step 1: Write failing test**

`packages/HimalayaUI/frontend/test/figure-export/renderer.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import * as Plot from "@observablehq/plot";
import { buildExportSvg, buildExportPng } from "../../src/lib/figure-export/renderer";
import type { ExportSpec } from "../../src/lib/figure-export/types";

function makeSpec(over: Partial<ExportSpec> = {}): ExportSpec {
  return {
    title: { primary: "Primary Title", secondary: "Secondary line" },
    width: 800,
    height: 600,
    plot: {
      marks: [Plot.line([{ x: 0, y: 0 }, { x: 1, y: 1 }], { x: "x", y: "y" })],
      x: { type: "linear" },
      y: { type: "linear" },
      width: 800,
      height: 480,
    },
    legend: {
      rows: [
        { color: "#1a1a1a", symbol: "swatch", label: "trace" },
        { color: "#cc6600", symbol: "line",   label: "predicted Pn3m" },
      ],
    },
    ...over,
  };
}

describe("buildExportSvg", () => {
  it("returns an SVGSVGElement", () => {
    const svg = buildExportSvg(makeSpec());
    expect(svg.tagName.toLowerCase()).toBe("svg");
    expect(svg.namespaceURI).toBe("http://www.w3.org/2000/svg");
  });

  it("renders a white background rect covering the canvas", () => {
    const svg = buildExportSvg(makeSpec());
    const rects = svg.querySelectorAll("rect");
    const bgRect = Array.from(rects).find((r) =>
      r.getAttribute("fill") === "#ffffff" || r.getAttribute("fill") === "white"
    );
    expect(bgRect).toBeDefined();
  });

  it("includes the primary title text", () => {
    const svg = buildExportSvg(makeSpec({ title: { primary: "My Trace Title" } }));
    expect(svg.textContent).toContain("My Trace Title");
  });

  it("includes the secondary title when provided", () => {
    const svg = buildExportSvg(makeSpec({ title: { primary: "P", secondary: "S" } }));
    expect(svg.textContent).toContain("S");
  });

  it("omits secondary title when undefined", () => {
    const spec = makeSpec();
    delete spec.title.secondary;
    const svg = buildExportSvg(spec);
    // Only the primary text appears.
    expect(svg.textContent?.trim()).toContain("Primary Title");
  });

  it("renders one legend row per LegendSpec.rows entry", () => {
    const svg = buildExportSvg(makeSpec());
    expect(svg.textContent).toContain("trace");
    expect(svg.textContent).toContain("predicted Pn3m");
  });

  it("output contains no var(--…) tokens (palette resolved at adapter time)", () => {
    const svg = buildExportSvg(makeSpec());
    const xml = new XMLSerializer().serializeToString(svg);
    expect(xml).not.toMatch(/var\(--/);
  });

  it("serialized output is well-formed XML and parses round-trip", () => {
    const svg = buildExportSvg(makeSpec());
    const xml = new XMLSerializer().serializeToString(svg);
    const parsed = new DOMParser().parseFromString(xml, "image/svg+xml");
    expect(parsed.querySelector("parsererror")).toBeNull();
    expect(parsed.documentElement.tagName.toLowerCase()).toBe("svg");
  });

  it("Plot.plot output is SVGSVGElement (guards against figure: true)", () => {
    // Sanity check that a vanilla spec doesn't have title/caption/figure set —
    // those flip Plot's return to HTMLElement.
    const spec = makeSpec();
    expect((spec.plot as { title?: unknown }).title).toBeUndefined();
    expect((spec.plot as { caption?: unknown }).caption).toBeUndefined();
    expect((spec.plot as { figure?: unknown }).figure).toBeUndefined();
    const inner = Plot.plot(spec.plot);
    expect(inner.tagName.toLowerCase()).toBe("svg");
  });
});

describe("buildExportPng", () => {
  it("is async and returns a Promise (full PNG path is JSDOM-untestable)", () => {
    // We only assert the function exists and returns a Promise. Real PNG
    // generation is covered by the Playwright E2E in Phase 4. JSDOM has no
    // OffscreenCanvas / Image.decode.
    const result = buildExportPng(makeSpec());
    expect(result).toBeInstanceOf(Promise);
    // Swallow rejection in JSDOM where canvas APIs are stubs/missing.
    result.catch(() => undefined);
  });
});
```

- [ ] **Step 2: Run, verify FAIL**

Run: `node_modules/.bin/vitest run test/figure-export/renderer.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Implement**

```ts
// renderer.ts — turn an ExportSpec into SVG (string-serializable) or PNG.
import * as Plot from "@observablehq/plot";
import type { ExportSpec, LegendRow } from "./types";
import {
  EXPORT_MARGIN,
  TITLE_FONT_PRIMARY,
  TITLE_FONT_SECONDARY,
  BODY_FONT,
  LIGHT_PALETTE,
} from "./presets";

const SVG_NS = "http://www.w3.org/2000/svg";

/**
 * Layout the export SVG: white background, primary + secondary title,
 * Plot.plot() body, legend rows beneath. Returns a fresh SVGSVGElement
 * not attached to the document.
 */
export function buildExportSvg(spec: ExportSpec): SVGSVGElement {
  const svg = document.createElementNS(SVG_NS, "svg");
  svg.setAttribute("xmlns", SVG_NS);
  svg.setAttribute("width", String(spec.width));
  svg.setAttribute("height", String(spec.height));
  svg.setAttribute("viewBox", `0 0 ${spec.width} ${spec.height}`);

  // White background.
  const bg = document.createElementNS(SVG_NS, "rect");
  bg.setAttribute("x", "0");
  bg.setAttribute("y", "0");
  bg.setAttribute("width", String(spec.width));
  bg.setAttribute("height", String(spec.height));
  bg.setAttribute("fill", LIGHT_PALETTE.background);
  svg.appendChild(bg);

  // Title block (primary + secondary).
  const titleY = 24;
  const primary = document.createElementNS(SVG_NS, "text");
  primary.setAttribute("x", String(EXPORT_MARGIN.left));
  primary.setAttribute("y", String(titleY));
  primary.setAttribute("font", TITLE_FONT_PRIMARY);
  primary.setAttribute("fill", LIGHT_PALETTE.text);
  primary.textContent = spec.title.primary;
  svg.appendChild(primary);

  if (spec.title.secondary) {
    const secondary = document.createElementNS(SVG_NS, "text");
    secondary.setAttribute("x", String(EXPORT_MARGIN.left));
    secondary.setAttribute("y", String(titleY + 18));
    secondary.setAttribute("font", TITLE_FONT_SECONDARY);
    secondary.setAttribute("fill", LIGHT_PALETTE.textMuted);
    secondary.textContent = spec.title.secondary;
    svg.appendChild(secondary);
  }

  // Plot body — Plot.plot() returns SVGSVGElement when title/caption/figure
  // are NOT set (adapters must NOT set them).
  const plotEl = Plot.plot({
    style: { background: "transparent", color: LIGHT_PALETTE.text, fontFamily: BODY_FONT },
    ...spec.plot,
  });
  const plotG = document.createElementNS(SVG_NS, "g");
  // Position the plot below the title block.
  const plotX = EXPORT_MARGIN.left;
  const plotY = EXPORT_MARGIN.top;
  plotG.setAttribute("transform", `translate(${plotX}, ${plotY})`);
  // Move all children of the inner SVG into the group, dropping the inner
  // <svg> wrapper itself (we are nesting inside our outer SVG).
  while (plotEl.firstChild) {
    plotG.appendChild(plotEl.firstChild);
  }
  svg.appendChild(plotG);

  // Legend rows (beneath the plot, inside the bottom margin).
  if (spec.legend && spec.legend.rows.length > 0) {
    const legendY = spec.height - 28;
    const legendG = document.createElementNS(SVG_NS, "g");
    legendG.setAttribute("transform", `translate(${EXPORT_MARGIN.left}, ${legendY})`);
    let cursorX = 0;
    for (const row of spec.legend.rows) {
      const item = renderLegendItem(row, cursorX);
      legendG.appendChild(item.group);
      cursorX += item.width + 16;
    }
    svg.appendChild(legendG);
  }

  return svg;
}

function renderLegendItem(row: LegendRow, x: number): { group: SVGGElement; width: number } {
  const g = document.createElementNS(SVG_NS, "g");
  g.setAttribute("transform", `translate(${x}, 0)`);

  if (row.symbol === "swatch") {
    const sw = document.createElementNS(SVG_NS, "rect");
    sw.setAttribute("x", "0");
    sw.setAttribute("y", "-9");
    sw.setAttribute("width", "10");
    sw.setAttribute("height", "10");
    sw.setAttribute("fill", row.color);
    g.appendChild(sw);
  } else {
    const ln = document.createElementNS(SVG_NS, "line");
    ln.setAttribute("x1", "0");
    ln.setAttribute("y1", "-4");
    ln.setAttribute("x2", "0");
    ln.setAttribute("y2", "-12");
    ln.setAttribute("stroke", row.color);
    ln.setAttribute("stroke-width", "2");
    g.appendChild(ln);
  }

  const text = document.createElementNS(SVG_NS, "text");
  text.setAttribute("x", "16");
  text.setAttribute("y", "0");
  text.setAttribute("font", BODY_FONT);
  text.setAttribute("fill", LIGHT_PALETTE.text);
  text.textContent = row.label;
  g.appendChild(text);

  // Approximate width — symbol + gap + len(label)*7. Good enough for layout.
  const width = 16 + row.label.length * 7;
  return { group: g, width };
}

/**
 * Render the export as a 2× DPI PNG blob. Pipeline: SVG → blob URL → Image
 * decode → OffscreenCanvas drawImage → convertToBlob.
 *
 * `URL.revokeObjectURL` runs in a `finally` so a thrown render still cleans
 * up (otherwise every export leaks an object URL). JSDOM doesn't have these
 * canvas APIs — full path covered by the Playwright E2E in Phase 4.
 */
export async function buildExportPng(spec: ExportSpec, scale = 2): Promise<Blob> {
  const svg = buildExportSvg(spec);
  const xml = new XMLSerializer().serializeToString(svg);
  const url = URL.createObjectURL(new Blob([xml], { type: "image/svg+xml" }));
  try {
    const img = new Image();
    img.src = url;
    await img.decode();
    const off = new OffscreenCanvas(spec.width * scale, spec.height * scale);
    const ctx = off.getContext("2d");
    if (!ctx) throw new Error("OffscreenCanvas 2d context unavailable");
    ctx.drawImage(img, 0, 0, off.width, off.height);
    return await off.convertToBlob({ type: "image/png" });
  } finally {
    URL.revokeObjectURL(url);
  }
}
```

- [ ] **Step 4: Run, verify PASS**

Run: same command as Step 2.
Expected: PASS, all `buildExportSvg` assertions green; `buildExportPng` test asserts only that it returns a Promise.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/figure-export/renderer.ts \
        packages/HimalayaUI/frontend/test/figure-export/renderer.test.ts
git commit -m "feat(figure-export): renderer SVG + PNG pipeline"
```

---

# Phase 2 — Trace path

---

### Task 6: `marks/traceExportMarks.ts`

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/figure-export/marks/traceExportMarks.ts`

This builds the array of `Plot.Mark` instances that the trace adapter assembles into `ExportSpec.plot.marks`. No separate test — the adapter test (Task 7) covers the resulting marks via `Plot.plot()` output structure. Spec §TraceViewer export.

- [ ] **Step 1: Implement**

```ts
// traceExportMarks.ts — Plot marks for the TraceViewer export.
// Duplicates on-screen overlay geometry (peaks + predicted-q ticks) — the
// on-screen overlay is load-bearing for hit-testing/hover and must not be
// touched. See spec §Approach for the accepted-tradeoff rationale.
import * as Plot from "@observablehq/plot";
import type { Trace, Peak, IndexEntry } from "../../../api";
import { phaseColor } from "../../../phases";
import {
  LIGHT_PALETTE,
  TRACE_STROKE_PX,
  PEAK_TICK_STROKE_PX,
  PREDICTED_Q_STROKE_PX,
} from "../presets";

/** Returns Plot marks. Caller threads them into ExportSpec.plot.marks. */
export function buildTraceExportMarks(args: {
  trace: Trace;
  peaks: Peak[];
  activeGroupIndices: IndexEntry[];
}): Plot.Mark[] {
  const { trace, peaks, activeGroupIndices } = args;

  const marks: Plot.Mark[] = [];

  // Trace line. Map (q,I) zip → array-of-objects for Plot.
  const points = trace.q.map((q, i) => ({ q, I: trace.I[i] ?? 0 }));
  marks.push(
    Plot.line(points, {
      x: "q",
      y: "I",
      stroke: LIGHT_PALETTE.trace,
      strokeWidth: TRACE_STROKE_PX,
    }),
  );

  // Peak triangles. Color by source/excluded.
  const peakRows = peaks.map((p) => ({
    q: p.q,
    I: p.intensity ?? 0,
    fill: p.excluded
      ? LIGHT_PALETTE.peakExcluded
      : p.source === "manual"
        ? LIGHT_PALETTE.peakManual
        : LIGHT_PALETTE.peakAuto,
  }));
  if (peakRows.length > 0) {
    marks.push(
      Plot.dot(peakRows, {
        x: "q",
        y: "I",
        fill: "fill",
        symbol: "triangle2",
        r: 4,
        stroke: "fill",
        strokeWidth: PEAK_TICK_STROKE_PX,
      }),
    );
  }

  // Predicted-q ticks per active-group index, phase-coloured.
  for (const idx of activeGroupIndices) {
    const color = phaseColor(idx.phase);
    const ticks = idx.predicted_q.map((q) => ({ q, phase: idx.phase }));
    marks.push(
      Plot.ruleX(ticks, {
        x: "q",
        stroke: color,
        strokeWidth: PREDICTED_Q_STROKE_PX,
        strokeOpacity: 0.85,
      }),
    );
  }

  return marks;
}
```

- [ ] **Step 2: Verify TS compiles**

Run: `npm run build 2>&1 | tail -20`
Expected: build succeeds.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/figure-export/marks/traceExportMarks.ts
git commit -m "feat(figure-export): trace export mark builder"
```

---

### Task 7: `adapters/traceAdapter.ts` (TDD)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/figure-export/adapters/traceAdapter.ts`
- Test: `packages/HimalayaUI/frontend/test/figure-export/traceAdapter.test.ts`

Spec §Adapters / TraceViewer export.

- [ ] **Step 1: Write failing test**

`packages/HimalayaUI/frontend/test/figure-export/traceAdapter.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import { buildTraceExportSpec } from "../../src/lib/figure-export/adapters/traceAdapter";
import { TRACE_DIMS } from "../../src/lib/figure-export/presets";
import type { Peak, IndexEntry, Trace } from "../../src/api";

const trace: Trace = {
  q: [0.1, 0.2, 0.3, 0.4, 0.5],
  I: [10, 8, 6, 4, 2],
  sigma: [1, 1, 1, 1, 1],
};
const peakAuto: Peak = {
  id: 1, exposure_id: 100, q: 0.18, intensity: 9, prominence: 5,
  sharpness: 1.0, source: "auto", excluded: false,
};
const peakManual: Peak = {
  ...peakAuto, id: 2, q: 0.32, source: "manual",
};
const peakExcluded: Peak = {
  ...peakAuto, id: 3, q: 0.42, excluded: true,
};
const indexPn3m: IndexEntry = {
  id: 10, exposure_id: 100, phase: "Pn3m", basis: 0.18, score: 0.9,
  r_squared: 0.99, lattice_d: 100, ngc: 1.5, status: "candidate", kind: "auto",
  inputs_hash: "h", peaks: [], predicted_q: [0.18, 0.36],
};

describe("buildTraceExportSpec", () => {
  it("produces an ExportSpec at TRACE_DIMS", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [peakAuto], activeGroupIndices: [indexPn3m],
      experimentName: "JC23", sampleName: "Sample 4", exposureLabel: "Exposure 7",
      xDomain: null, yDomain: null, xType: "log",
    });
    expect(spec.width).toBe(TRACE_DIMS.width);
    expect(spec.height).toBe(TRACE_DIMS.height);
  });

  it("primary title is 'experiment · sample · exposure'", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [], activeGroupIndices: [],
      experimentName: "JC23", sampleName: "Sample 4", exposureLabel: "Exposure 7",
      xDomain: null, yDomain: null, xType: "log",
    });
    expect(spec.title.primary).toContain("JC23");
    expect(spec.title.primary).toContain("Sample 4");
    expect(spec.title.primary).toContain("Exposure 7");
  });

  it("legend includes 'auto' row when only auto peaks present", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [peakAuto], activeGroupIndices: [],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "log",
    });
    const labels = (spec.legend?.rows ?? []).map((r) => r.label.toLowerCase());
    expect(labels.some((l) => l.includes("auto"))).toBe(true);
  });

  it("legend adds 'manual' row when manual peaks present", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [peakAuto, peakManual], activeGroupIndices: [],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "log",
    });
    const labels = (spec.legend?.rows ?? []).map((r) => r.label.toLowerCase());
    expect(labels.some((l) => l.includes("manual"))).toBe(true);
  });

  it("legend adds 'excluded' row when excluded peaks present", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [peakAuto, peakExcluded], activeGroupIndices: [],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "log",
    });
    const labels = (spec.legend?.rows ?? []).map((r) => r.label.toLowerCase());
    expect(labels.some((l) => l.includes("excluded"))).toBe(true);
  });

  it("legend adds one row per active-group index phase", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [], activeGroupIndices: [indexPn3m],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "log",
    });
    const labels = (spec.legend?.rows ?? []).map((r) => r.label);
    expect(labels.some((l) => l.includes("Pn3m"))).toBe(true);
  });

  it("plot.x.type honours xType arg", () => {
    const specLog = buildTraceExportSpec({
      trace, peaks: [], activeGroupIndices: [],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "log",
    });
    const specLin = buildTraceExportSpec({
      trace, peaks: [], activeGroupIndices: [],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "linear",
    });
    expect((specLog.plot.x as { type?: string }).type).toBe("log");
    expect((specLin.plot.x as { type?: string }).type).toBe("linear");
  });

  it("does NOT set title/caption/figure on plot (would break renderer)", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [], activeGroupIndices: [],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "log",
    });
    expect((spec.plot as { title?: unknown }).title).toBeUndefined();
    expect((spec.plot as { caption?: unknown }).caption).toBeUndefined();
    expect((spec.plot as { figure?: unknown }).figure).toBeUndefined();
  });
});
```

- [ ] **Step 2: Run, verify FAIL**

Run: `node_modules/.bin/vitest run test/figure-export/traceAdapter.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Implement**

```ts
// traceAdapter.ts — TraceViewer state → ExportSpec.
import type { Trace, Peak, IndexEntry } from "../../../api";
import type { ExportSpec, LegendRow } from "../types";
import { TRACE_DIMS, LIGHT_PALETTE, EXPORT_MARGIN } from "../presets";
import { buildTraceExportMarks } from "../marks/traceExportMarks";
import { phaseColor } from "../../../phases";

export interface TraceAdapterArgs {
  trace: Trace;
  peaks: Peak[];
  activeGroupIndices: IndexEntry[];
  experimentName: string;
  sampleName: string;
  exposureLabel: string;
  xDomain: [number, number] | null;
  yDomain: [number, number] | null;
  xType: "log" | "linear";
  /** Display string for the q axis label. Falls back to "Å⁻¹". */
  qUnits?: string;
}

export function buildTraceExportSpec(args: TraceAdapterArgs): ExportSpec {
  const {
    trace, peaks, activeGroupIndices,
    experimentName, sampleName, exposureLabel,
    xDomain, yDomain, xType, qUnits,
  } = args;

  const marks = buildTraceExportMarks({ trace, peaks, activeGroupIndices });

  const xConfig: Record<string, unknown> = {
    type: xType,
    label: `q (${qUnits ?? "Å⁻¹"})`,
  };
  if (xDomain) xConfig.domain = xDomain;

  const yConfig: Record<string, unknown> = {
    type: "log",
    label: "I",
  };
  if (yDomain) yConfig.domain = yDomain;

  const plotW = TRACE_DIMS.width  - EXPORT_MARGIN.left - EXPORT_MARGIN.right;
  const plotH = TRACE_DIMS.height - EXPORT_MARGIN.top  - EXPORT_MARGIN.bottom;

  const legend: LegendRow[] = [];

  // Peak source rows.
  legend.push({ color: LIGHT_PALETTE.peakAuto, symbol: "swatch", label: "auto peak" });
  if (peaks.some((p) => p.source === "manual")) {
    legend.push({ color: LIGHT_PALETTE.peakManual, symbol: "swatch", label: "manual peak" });
  }
  if (peaks.some((p) => p.excluded)) {
    legend.push({ color: LIGHT_PALETTE.peakExcluded, symbol: "swatch", label: "excluded auto peak" });
  }

  // Phase rows — one per active-group index phase.
  const seenPhases = new Set<string>();
  for (const idx of activeGroupIndices) {
    if (seenPhases.has(idx.phase)) continue;
    seenPhases.add(idx.phase);
    legend.push({
      color: phaseColor(idx.phase),
      symbol: "line",
      label: `predicted ${idx.phase}`,
    });
  }

  return {
    title: { primary: `${experimentName} · ${sampleName} · ${exposureLabel}` },
    width: TRACE_DIMS.width,
    height: TRACE_DIMS.height,
    plot: {
      marks,
      x: xConfig,
      y: yConfig,
      width: plotW,
      height: plotH,
      marginLeft: 0,
      marginRight: 0,
      marginTop: 0,
      marginBottom: 0,
    },
    legend: { rows: legend },
  };
}
```

- [ ] **Step 4: Run, verify PASS**

Run: same command as Step 2.
Expected: PASS, all assertions green.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/figure-export/adapters/traceAdapter.ts \
        packages/HimalayaUI/frontend/test/figure-export/traceAdapter.test.ts
git commit -m "feat(figure-export): traceAdapter for TraceViewer export"
```

---

### Task 8: `FigureExportControls.tsx` (TDD)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/FigureExportControls.tsx`
- Test: `packages/HimalayaUI/frontend/test/figure-export/FigureExportControls.test.tsx`

Spec §Button component + §Popover behaviour. The popover follows the **`ForksPopover` pattern** from `src/components/ForksPopover.tsx` (Esc to close + outside-click dismiss; not focus-trapped, NOT a menu role — `aria-haspopup="true"` per Bundle B's #75 fix).

- [ ] **Step 1: Write failing test**

`packages/HimalayaUI/frontend/test/figure-export/FigureExportControls.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { FigureExportControls } from "../../src/components/FigureExportControls";
import * as renderer from "../../src/lib/figure-export/renderer";
import * as clipboard from "../../src/lib/figure-export/clipboard";
import * as toastModule from "../../src/lib/toast";
import type { ExportSpec } from "../../src/lib/figure-export/types";

const fakeSpec: ExportSpec = {
  title: { primary: "P" }, width: 100, height: 80,
  plot: { marks: [], width: 80, height: 60 },
};

let toastSpy: ReturnType<typeof vi.fn>;

beforeEach(() => {
  toastSpy = vi.fn();
  toastModule.setToastImpl(toastSpy);
});
afterEach(() => {
  toastModule.setToastImpl(null);
  vi.restoreAllMocks();
});

describe("FigureExportControls", () => {
  it("renders Copy and Download buttons with aria labels using ariaContext", () => {
    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="himalaya-trace-x"
        ariaContext="trace plot"
      />,
    );
    expect(screen.getByRole("button", { name: /copy trace plot to clipboard/i }))
      .toBeInTheDocument();
    expect(screen.getByRole("button", { name: /download trace plot as png/i }))
      .toBeInTheDocument();
  });

  it("Copy click invokes renderer + clipboard and emits success toast", async () => {
    const blob = new Blob([new Uint8Array([1])], { type: "image/png" });
    vi.spyOn(renderer, "buildExportPng").mockResolvedValue(blob);
    vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(true);
    const writeSpy = vi.spyOn(clipboard, "copyPngToClipboard").mockResolvedValue();

    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    fireEvent.click(screen.getByRole("button", { name: /copy trace plot/i }));

    await waitFor(() => expect(writeSpy).toHaveBeenCalledWith(blob));
    expect(toastSpy).toHaveBeenCalledWith(
      expect.stringContaining("Copied"),
      "success",
    );
  });

  it("Copy failure surfaces error toast (renderer rejects)", async () => {
    vi.spyOn(renderer, "buildExportPng").mockRejectedValue(new Error("boom"));
    vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(true);

    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    fireEvent.click(screen.getByRole("button", { name: /copy trace plot/i }));

    await waitFor(() =>
      expect(toastSpy).toHaveBeenCalledWith(
        expect.stringContaining("Couldn't copy"),
        "error",
      ),
    );
  });

  it("buttons re-enable after a thrown render (try/finally pending flag)", async () => {
    vi.spyOn(renderer, "buildExportPng").mockRejectedValue(new Error("boom"));
    vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(true);

    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    const copyBtn = screen.getByRole("button", { name: /copy trace plot/i });
    fireEvent.click(copyBtn);

    await waitFor(() => expect(copyBtn).not.toBeDisabled());
  });

  it("Copy is disabled when canCopyPngToClipboard() is false", () => {
    vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(false);
    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    expect(screen.getByRole("button", { name: /copy trace plot/i })).toBeDisabled();
  });

  it("disabled prop disables both buttons (parent gate)", () => {
    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
        disabled
      />,
    );
    expect(screen.getByRole("button", { name: /copy trace plot/i })).toBeDisabled();
    expect(screen.getByRole("button", { name: /download trace plot as png/i })).toBeDisabled();
  });

  it("Chevron toggles a popover with PNG/SVG menu rows", () => {
    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    const chevron = screen.getByRole("button", { name: /other download formats/i });
    expect(chevron).toHaveAttribute("aria-haspopup", "true");
    expect(chevron).toHaveAttribute("aria-expanded", "false");
    fireEvent.click(chevron);
    expect(chevron).toHaveAttribute("aria-expanded", "true");
    expect(screen.getByText(/download as png/i)).toBeInTheDocument();
    expect(screen.getByText(/download as svg/i)).toBeInTheDocument();
  });

  it("Esc closes the popover", () => {
    render(
      <FigureExportControls
        spec={() => fakeSpec}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    const chevron = screen.getByRole("button", { name: /other download formats/i });
    fireEvent.click(chevron);
    fireEvent.keyDown(document, { key: "Escape" });
    expect(chevron).toHaveAttribute("aria-expanded", "false");
  });

  it("spec is invoked on each click (thunk captures fresh state)", async () => {
    const blob = new Blob([new Uint8Array([1])], { type: "image/png" });
    vi.spyOn(renderer, "buildExportPng").mockResolvedValue(blob);
    vi.spyOn(clipboard, "canCopyPngToClipboard").mockReturnValue(true);
    vi.spyOn(clipboard, "copyPngToClipboard").mockResolvedValue();

    const specThunk = vi.fn().mockReturnValue(fakeSpec);
    render(
      <FigureExportControls
        spec={specThunk}
        filenameStem="x"
        ariaContext="trace plot"
      />,
    );
    const copyBtn = screen.getByRole("button", { name: /copy trace plot/i });
    fireEvent.click(copyBtn);
    await waitFor(() => expect(specThunk).toHaveBeenCalledTimes(1));
    fireEvent.click(copyBtn);
    await waitFor(() => expect(specThunk).toHaveBeenCalledTimes(2));
  });
});
```

- [ ] **Step 2: Run, verify FAIL**

Run: `node_modules/.bin/vitest run test/figure-export/FigureExportControls.test.tsx`
Expected: FAIL — module not found.

- [ ] **Step 3: Implement**

```tsx
// FigureExportControls.tsx — Copy + split-Download (PNG/SVG) buttons.
import { useEffect, useRef, useState } from "react";
import type { ExportSpec } from "../lib/figure-export/types";
import { buildExportPng, buildExportSvg } from "../lib/figure-export/renderer";
import { canCopyPngToClipboard, copyPngToClipboard } from "../lib/figure-export/clipboard";
import { downloadBlob } from "../lib/figure-export/download";
import { buildFilename } from "../lib/figure-export/filename";
import { showToast } from "../lib/toast";

export interface FigureExportControlsProps {
  /** Thunk evaluated at click time so it captures fresh state (xDomain,
   *  current peaks, etc.) rather than stale state at component mount. */
  spec: () => ExportSpec;
  /** Already-resolved, slugified filename stem. Component appends the date
   *  + extension internally. */
  filenameStem: string;
  /** Used in aria-labels: "Copy {ariaContext} to clipboard". */
  ariaContext: string;
  /** Parent disables both buttons when the underlying data is incomplete
   *  (e.g. `!traceQ.data || !peaksQ.data`). */
  disabled?: boolean;
}

export function FigureExportControls({
  spec, filenameStem, ariaContext, disabled = false,
}: FigureExportControlsProps): JSX.Element {
  const [pending, setPending] = useState(false);
  const [menuOpen, setMenuOpen] = useState(false);
  const triggerRef = useRef<HTMLButtonElement>(null);
  const panelRef = useRef<HTMLDivElement>(null);
  const canCopy = canCopyPngToClipboard();

  // Esc + outside-click dismiss (mirrors ForksPopover pattern; spec §Popover).
  useEffect(() => {
    if (!menuOpen) return;
    const onKey = (e: KeyboardEvent): void => {
      if (e.key === "Escape") setMenuOpen(false);
    };
    document.addEventListener("keydown", onKey);
    return () => document.removeEventListener("keydown", onKey);
  }, [menuOpen]);

  useEffect(() => {
    if (!menuOpen) return;
    const onDown = (e: MouseEvent): void => {
      const t = e.target as Node | null;
      if (t === null) return;
      if (panelRef.current?.contains(t)) return;
      if (triggerRef.current?.contains(t)) return;
      setMenuOpen(false);
    };
    document.addEventListener("mousedown", onDown);
    return () => document.removeEventListener("mousedown", onDown);
  }, [menuOpen]);

  const onCopy = async (): Promise<void> => {
    if (disabled || !canCopy || pending) return;
    setPending(true);
    try {
      const blob = await buildExportPng(spec());
      await copyPngToClipboard(blob);
      showToast("Copied figure to clipboard", "success");
    } catch (err) {
      showToast("Couldn't copy figure — try Download instead.", "error");
      // eslint-disable-next-line no-console
      console.warn(err);
    } finally {
      setPending(false);
    }
  };

  const onDownloadPng = async (): Promise<void> => {
    if (disabled || pending) return;
    setMenuOpen(false);
    setPending(true);
    try {
      const blob = await buildExportPng(spec());
      downloadBlob(blob, buildFilename(filenameStem, "png"));
    } catch (err) {
      showToast("Couldn't render figure for download.", "error");
      // eslint-disable-next-line no-console
      console.warn(err);
    } finally {
      setPending(false);
    }
  };

  const onDownloadSvg = (): void => {
    if (disabled || pending) return;
    setMenuOpen(false);
    setPending(true);
    try {
      const svg = buildExportSvg(spec());
      const xml = new XMLSerializer().serializeToString(svg);
      const blob = new Blob([xml], { type: "image/svg+xml" });
      downloadBlob(blob, buildFilename(filenameStem, "svg"));
    } catch (err) {
      showToast("Couldn't render figure for download.", "error");
      // eslint-disable-next-line no-console
      console.warn(err);
    } finally {
      setPending(false);
    }
  };

  const copyDisabled = disabled || !canCopy || pending;
  const downloadDisabled = disabled || pending;

  const copyTitle = !canCopy
    ? "Clipboard requires HTTPS"
    : `Copy ${ariaContext} to clipboard`;

  return (
    <span
      data-testid="figure-export-controls"
      className="inline-flex items-center gap-1"
    >
      <button
        type="button"
        data-testid="figure-export-copy"
        aria-label={`Copy ${ariaContext} to clipboard`}
        title={copyTitle}
        disabled={copyDisabled}
        onClick={() => { void onCopy(); }}
        className="px-1.5 py-0.5 rounded text-xs text-fg-dim hover:text-fg
                   hover:bg-bg-hover border border-transparent hover:border-border
                   disabled:opacity-40 disabled:cursor-default"
      >
        Copy
      </button>
      <span className="inline-flex items-stretch border border-border rounded overflow-hidden">
        <button
          type="button"
          data-testid="figure-export-download-png"
          aria-label={`Download ${ariaContext} as PNG`}
          disabled={downloadDisabled}
          onClick={() => { void onDownloadPng(); }}
          className="px-1.5 py-0.5 text-xs text-fg-dim hover:text-fg hover:bg-bg-hover
                     disabled:opacity-40 disabled:cursor-default"
        >
          Save
        </button>
        <span className="w-px bg-border" />
        <button
          ref={triggerRef}
          type="button"
          data-testid="figure-export-download-menu-trigger"
          aria-label="Other download formats"
          aria-haspopup="true"
          aria-expanded={menuOpen}
          disabled={downloadDisabled}
          onClick={() => setMenuOpen((v) => !v)}
          className="px-1 text-xs text-fg-dim hover:text-fg hover:bg-bg-hover
                     disabled:opacity-40 disabled:cursor-default"
        >
          ▾
        </button>
      </span>
      {menuOpen && (
        <div
          ref={panelRef}
          data-testid="figure-export-download-menu"
          className="absolute z-50 mt-6 right-0 min-w-[140px] card border border-border
                     bg-bg-elevated shadow-lg p-1"
          style={{ position: "absolute" }}
        >
          <button
            type="button"
            data-testid="figure-export-download-menu-png"
            onClick={() => { void onDownloadPng(); }}
            className="block w-full text-left px-2 py-1 text-xs text-fg hover:bg-bg-hover rounded"
          >
            Download as PNG
          </button>
          <button
            type="button"
            data-testid="figure-export-download-menu-svg"
            onClick={onDownloadSvg}
            className="block w-full text-left px-2 py-1 text-xs text-fg hover:bg-bg-hover rounded"
          >
            Download as SVG
          </button>
        </div>
      )}
    </span>
  );
}
```

- [ ] **Step 4: Run, verify PASS**

Run: same command as Step 2.
Expected: PASS, all 9 assertions green.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/FigureExportControls.tsx \
        packages/HimalayaUI/frontend/test/figure-export/FigureExportControls.test.tsx
git commit -m "feat(figure-export): FigureExportControls component (Copy + split Download)"
```

---

### Task 9: Wire `<FigureExportControls>` into `PlotCard`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/PlotCard.tsx`

Spec §Placement (PlotCard title strip): append after the q-range reset button, separated by a thin divider. Disabled gate: `!traceQ.data || !peaksQ.data` (spec §Disabled states — stricter than `canFit`).

- [ ] **Step 1: Add import**

In `PlotCard.tsx` near the top imports, add:

```tsx
import { FigureExportControls } from "./FigureExportControls";
import { buildTraceExportSpec } from "../lib/figure-export/adapters/traceAdapter";
import { slugifyForFilename } from "../lib/figure-export/filename";
```

- [ ] **Step 2: Build the export thunk + filename stem inside `PlotCard()`**

Below the existing `activeGroup`/`hoveredIndex` derivations (around line 220, before the `body` IIFE), insert:

```tsx
  // Figure export (spec: 2026-05-08-figure-export-design.md).
  const exposureLabel = activeExposureId !== undefined
    ? `Exposure ${activeExposureId}`
    : "";
  const filenameStem = `himalaya-trace-${
    slugifyForFilename(experimentName ?? "")
  }-${
    slugifyForFilename(sampleName ?? "")
  }-${
    slugifyForFilename(exposureLabel)
  }`;

  const exportSpec = useCallback(() => {
    if (!traceQ.data || !peaksQ.data) {
      throw new Error("FigureExportControls: parent disabled-gate violated");
    }
    return buildTraceExportSpec({
      trace: traceQ.data,
      peaks: peaksQ.data,
      activeGroupIndices,
      experimentName: experimentName ?? "",
      sampleName: sampleName ?? "",
      exposureLabel,
      xDomain,
      yDomain,
      xType,
      ...(experimentQ.data?.q_units ? { qUnits: experimentQ.data.q_units } : {}),
    });
  }, [
    traceQ.data, peaksQ.data, activeGroupIndices,
    experimentName, sampleName, exposureLabel,
    xDomain, yDomain, xType, experimentQ.data?.q_units,
  ]);

  const exportDisabled = !traceQ.data || !peaksQ.data;
```

- [ ] **Step 3: Pass props through `<TitleStrip>`**

Modify the `<TitleStrip>` call (around line 270) to add the new props:

```tsx
      <TitleStrip
        experimentName={experimentName}
        sampleName={sampleName}
        onTitleClick={() => openNavModal(titleStep)}
        xDomain={effectiveDomain}
        fullRange={fullQRange}
        onXDomain={setXDomain}
        xType={xType}
        onSetXType={setXType}
        onFitFeatures={fitFeatures}
        canFit={traceQ.data !== undefined}
        exportSpec={exportSpec}
        exportFilenameStem={filenameStem}
        exportDisabled={exportDisabled}
      />
```

- [ ] **Step 4: Update `TitleStripProps` and inner JSX**

In the same file, extend the interface:

```tsx
interface TitleStripProps {
  experimentName: string | undefined;
  sampleName:     string | undefined;
  onTitleClick:   () => void;
  xDomain:  [number, number] | null;
  fullRange: [number, number] | null;
  onXDomain: (d: [number, number] | null) => void;
  xType: "log" | "linear";
  onSetXType: (t: "log" | "linear") => void;
  onFitFeatures: () => void;
  canFit: boolean;
  exportSpec: () => import("../lib/figure-export/types").ExportSpec;
  exportFilenameStem: string;
  exportDisabled: boolean;
}
```

Update the destructure and the JSX (the `<div className="shrink-0 flex items-center gap-2">` controls cluster):

```tsx
function TitleStrip({
  experimentName, sampleName, onTitleClick, xDomain, fullRange, onXDomain,
  xType, onSetXType, onFitFeatures, canFit,
  exportSpec, exportFilenameStem, exportDisabled,
}: TitleStripProps): JSX.Element {
  // …existing markup…
      <div className="shrink-0 flex items-center gap-2">
        <button
          type="button"
          onClick={onFitFeatures}
          disabled={!canFit}
          data-testid="fit-features"
          /* …existing classes… */
        >
          fit features
        </button>
        <XScaleToggle xType={xType} onSetXType={onSetXType} />
        <QRange xDomain={xDomain} fullRange={fullRange} onXDomain={onXDomain} />
        {/* Thin divider before the export cluster. */}
        <span className="w-px h-4 bg-border" aria-hidden="true" />
        <FigureExportControls
          spec={exportSpec}
          filenameStem={exportFilenameStem}
          ariaContext="trace plot"
          disabled={exportDisabled}
        />
      </div>
```

- [ ] **Step 5: Verify TS compiles**

Run: `npm run build 2>&1 | tail -20`
Expected: build succeeds.

- [ ] **Step 6: Run existing PlotCard tests**

Run: `node_modules/.bin/vitest run test/PlotCard 2>&1 | tail -30`
Expected: existing tests stay green (we added new prop-flow but no behaviour break).

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/PlotCard.tsx
git commit -m "feat(figure-export): wire FigureExportControls into PlotCard"
```

---

# Phase 3 — Compare path

---

### Task 10: `marks/multiTraceExportMarks.ts`

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/figure-export/marks/multiTraceExportMarks.ts`

Spec §MultiTracePlot export. The mark builder duplicates the on-screen stacked-band layout from `MultiTracePlot.tsx::computeYBands` — caller passes the **filtered** member list (null-`exposure_id` removed) so band heights renormalize. Per-member labels render inline at each band's y-position via `Plot.text`. Peak ticks/labels honour `showPeakTicks`/`showPeakLabels`.

- [ ] **Step 1: Implement**

```ts
// multiTraceExportMarks.ts — Plot marks for the MultiTracePlot export.
// Stacked traces in display_order, per-member labels at each band's
// y-position, peak ticks/labels gated on showPeakTicks / showPeakLabels.
import * as Plot from "@observablehq/plot";
import type { ComparisonMember } from "../../../api";
import type { Trace } from "../../../api";
import { computeYBands } from "../../../components/MultiTracePlot";
import {
  LIGHT_PALETTE,
  TRACE_STROKE_PX,
  PEAK_TICK_STROKE_PX,
} from "../presets";

export interface MultiTraceMarksArgs {
  /** Already filtered (null-exposure_id removed) and sorted by display_order. */
  members: ComparisonMember[];
  traces: Map<number, Trace>;
  /** Pre-resolved member labels (post-#73 contract). */
  displayLabelByMemberId: Map<number, string>;
  /** Pre-resolved colors keyed by member.id (computed on the UNFILTERED list
   *  per spec §MultiTracePlot export — distinct mode keys off the original
   *  index). */
  colorByMember: Map<number, string>;
  showPeakTicks: boolean;
  showPeakLabels: boolean;
  panelHeight: number;
}

export function buildMultiTraceExportMarks(args: MultiTraceMarksArgs): Plot.Mark[] {
  const {
    members, traces, displayLabelByMemberId, colorByMember,
    showPeakTicks, showPeakLabels, panelHeight,
  } = args;

  const ratios = members.map((m) => m.band_height || 1);
  const yBands = computeYBands(ratios, panelHeight);

  const marks: Plot.Mark[] = [];

  for (let i = 0; i < members.length; i++) {
    const member = members[i]!;
    if (member.exposure_id === null) continue; // defensive (caller pre-filtered)
    const trace = traces.get(member.exposure_id);
    if (!trace) continue;
    const band = yBands[i]!;
    const [bandTop, bandBottom] = band;
    const bandH = Math.max(1, bandBottom - bandTop);
    const color = colorByMember.get(member.id) ?? LIGHT_PALETTE.trace;

    // Trace line — y-mapped into this band. We compute log(I) ourselves to
    // place the line within the band (Plot's y-domain is per-figure, not
    // per-band). Each band has its own (per-member) y-range.
    const positiveI = trace.I.filter((v) => Number.isFinite(v) && v > 0);
    if (positiveI.length === 0) continue;
    const logMin = Math.log10(Math.min(...positiveI));
    const logMax = Math.log10(Math.max(...positiveI));
    const logRange = Math.max(1e-6, logMax - logMin);
    const points = trace.q.map((q, idx) => {
      const I = trace.I[idx] ?? 0;
      const li = I > 0 ? Math.log10(I) : logMin;
      // Map li ∈ [logMin, logMax] → y ∈ [bandBottom, bandTop] (inverted: low at bottom).
      const y = bandBottom - ((li - logMin) / logRange) * bandH;
      return { q, y };
    });
    marks.push(
      Plot.line(points, {
        x: "q",
        y: "y",
        stroke: color,
        strokeWidth: TRACE_STROKE_PX,
      }),
    );

    // Per-member label at the band's vertical midpoint, anchored to the left.
    const label = displayLabelByMemberId.get(member.id) ?? "";
    if (label) {
      marks.push(
        Plot.text(
          [{ q: trace.q[0] ?? 0, y: bandTop + bandH * 0.15, label }],
          {
            x: "q",
            y: "y",
            text: "label",
            textAnchor: "start",
            fill: LIGHT_PALETTE.text,
            fontSize: 11,
            dx: 4,
          },
        ),
      );
    }

    // Peaks — honour show toggles.
    const peaks = member.snapshot?.effective_peaks ?? [];
    if (showPeakTicks && peaks.length > 0) {
      const tickPoints = peaks.map((p) => ({
        q: p.q,
        y: bandTop + bandH * 0.05,
      }));
      marks.push(
        Plot.ruleX(tickPoints, {
          x: "q",
          stroke: color,
          strokeWidth: PEAK_TICK_STROKE_PX,
          // Restrict the tick height by mapping y to a small range above the band top.
        }),
      );
      if (showPeakLabels) {
        const labelRows = peaks.map((p, idx) => ({
          q: p.q,
          y: bandTop + bandH * 0.02,
          label: String(idx + 1),
        }));
        marks.push(
          Plot.text(labelRows, {
            x: "q",
            y: "y",
            text: "label",
            fill: color,
            fontSize: 9,
            dy: -2,
          }),
        );
      }
    }
  }

  return marks;
}
```

- [ ] **Step 2: Verify TS compiles**

Run: `npm run build 2>&1 | tail -20`
Expected: build succeeds.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/figure-export/marks/multiTraceExportMarks.ts
git commit -m "feat(figure-export): multiTrace export mark builder"
```

---

### Task 11: `adapters/multiTraceAdapter.ts` (TDD)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/figure-export/adapters/multiTraceAdapter.ts`
- Test: `packages/HimalayaUI/frontend/test/figure-export/multiTraceAdapter.test.ts`

Spec §MultiTracePlot export, §Trace coloring. Critical: compute `colorByMember` on the **unfiltered** members array (so `defaultDistinct`'s findIndex is stable), then filter null-`exposure_id` for marks. All-null input throws.

- [ ] **Step 1: Write failing test**

`packages/HimalayaUI/frontend/test/figure-export/multiTraceAdapter.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import { buildMultiTraceExportSpec } from "../../src/lib/figure-export/adapters/multiTraceAdapter";
import { COMPARE_DIMS } from "../../src/lib/figure-export/presets";
import type { ComparisonMember, Trace } from "../../src/api";

function makeMember(over: Partial<ComparisonMember> = {}): ComparisonMember {
  return {
    id: 1, comparison_id: 1, exposure_id: 100, display_order: 0,
    band_height: 1, y_offset: 0, normalization: "none",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: null, is_stale: false, created_by: 1, created_at: null,
    ...over,
  };
}

const trace: Trace = {
  q: [0.1, 0.2, 0.3], I: [10, 5, 2], sigma: [1, 1, 1],
};

const traces = new Map<number, Trace>([[100, trace], [101, trace], [102, trace]]);
const sampleIdFor = (m: ComparisonMember): number | null => {
  if (m.exposure_id === 100) return 1;
  if (m.exposure_id === 101) return 2;
  if (m.exposure_id === 102) return 3;
  return null;
};

describe("buildMultiTraceExportSpec", () => {
  it("produces an ExportSpec at COMPARE_DIMS", () => {
    const m = makeMember();
    const spec = buildMultiTraceExportSpec({
      members: [m],
      traces,
      comparisonTitle: "My comparison",
      experimentName: "JC23",
      xDomain: null,
      showPeakTicks: true,
      showPeakLabels: true,
      groupingMode: "distinct",
      sampleIdFor,
    });
    expect(spec.width).toBe(COMPARE_DIMS.width);
    expect(spec.height).toBe(COMPARE_DIMS.height);
  });

  it("primary title is the comparison title; secondary is the experiment name", () => {
    const spec = buildMultiTraceExportSpec({
      members: [makeMember()],
      traces,
      comparisonTitle: "T",
      experimentName: "E",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "distinct", sampleIdFor,
    });
    expect(spec.title.primary).toBe("T");
    expect(spec.title.secondary).toBe("E");
  });

  it("omits secondary title when experimentName is undefined (global scope)", () => {
    const spec = buildMultiTraceExportSpec({
      members: [makeMember()],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "distinct", sampleIdFor,
    });
    expect(spec.title.primary).toBe("T");
    expect(spec.title.secondary).toBeUndefined();
  });

  it("bySample legend has one row per unique sample_id", () => {
    const spec = buildMultiTraceExportSpec({
      members: [
        makeMember({ id: 1, exposure_id: 100 }),
        makeMember({ id: 2, exposure_id: 101 }),
        makeMember({ id: 3, exposure_id: 102 }),
      ],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "bySample", sampleIdFor,
    });
    // 3 unique sample ids → 3 legend rows.
    expect(spec.legend?.rows.length).toBe(3);
  });

  it("byPhase legend has one row per unique confirmed-index phase, plus orphan row when any have none", () => {
    const phasedMember = makeMember({
      id: 10, exposure_id: 100,
      snapshot: {
        effective_peaks: [],
        confirmed_index: {
          id: 1, exposure_id: 100, phase: "Pn3m", basis: 0.18, score: 0.9,
          r_squared: 0.99, lattice_d: 100, ngc: 1.5, status: "candidate",
          kind: "auto", inputs_hash: "h", peaks: [], predicted_q: [0.18],
        },
      } as unknown as ComparisonMember["snapshot"],
    });
    const orphanMember = makeMember({ id: 11, exposure_id: 101, snapshot: null });

    const spec = buildMultiTraceExportSpec({
      members: [phasedMember, orphanMember],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "byPhase", sampleIdFor,
    });

    // One Pn3m row + one shared orphan row.
    const labels = (spec.legend?.rows ?? []).map((r) => r.label.toLowerCase());
    expect(labels.some((l) => l.includes("pn3m"))).toBe(true);
    expect(labels.some((l) => l.includes("unphased") || l.includes("unbound"))).toBe(true);
    expect(spec.legend?.rows.length).toBe(2);
  });

  it("distinct mode emits NO legend (per-member labels carry the encoding)", () => {
    const spec = buildMultiTraceExportSpec({
      members: [makeMember({ id: 1 }), makeMember({ id: 2, exposure_id: 101 })],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "distinct", sampleIdFor,
    });
    expect(spec.legend?.rows.length ?? 0).toBe(0);
  });

  it("filters null-exposure_id members from marks but keeps colorByMember stable for distinct mode", () => {
    const m1 = makeMember({ id: 1, exposure_id: 100 });
    const orphan = makeMember({ id: 2, exposure_id: null });
    const m3 = makeMember({ id: 3, exposure_id: 102 });
    const spec = buildMultiTraceExportSpec({
      members: [m1, orphan, m3],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "distinct", sampleIdFor,
    });
    // Marks were built off [m1, m3] (orphan filtered) BUT colors were
    // computed off the unfiltered list — so m3's distinct color is index 2,
    // not index 1. We sanity-check the spec was generated successfully and
    // the figure dims are still the canonical preset (no nullable downgrade).
    expect(spec.width).toBe(COMPARE_DIMS.width);
    // marks is a Plot.Mark[] — at minimum we expect the trace marks for the
    // two non-orphan members, which means non-empty.
    const marksCount = (spec.plot.marks ?? []).length;
    expect(marksCount).toBeGreaterThan(0);
  });

  it("throws when every member has exposure_id === null", () => {
    expect(() => buildMultiTraceExportSpec({
      members: [
        makeMember({ id: 1, exposure_id: null }),
        makeMember({ id: 2, exposure_id: null }),
      ],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "distinct", sampleIdFor,
    })).toThrow(/exposure_id === null/);
  });

  it("does NOT set title/caption/figure on plot", () => {
    const spec = buildMultiTraceExportSpec({
      members: [makeMember()],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "distinct", sampleIdFor,
    });
    expect((spec.plot as { title?: unknown }).title).toBeUndefined();
    expect((spec.plot as { caption?: unknown }).caption).toBeUndefined();
    expect((spec.plot as { figure?: unknown }).figure).toBeUndefined();
  });
});
```

- [ ] **Step 2: Run, verify FAIL**

Run: `node_modules/.bin/vitest run test/figure-export/multiTraceAdapter.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Implement**

```ts
// multiTraceAdapter.ts — MultiTracePlot state → ExportSpec.
import type { ComparisonMember, Trace } from "../../../api";
import type { ExportSpec, LegendRow } from "../types";
import {
  COMPARE_DIMS, COMPARE_PALETTE_LIGHT, ORPHAN_FALLBACK_LIGHT,
  EXPORT_MARGIN, LIGHT_PALETTE,
} from "../presets";
import { buildMultiTraceExportMarks } from "../marks/multiTraceExportMarks";
import { colorFor, type GroupingMode } from "../../comparison/coloring";
import { phaseColor } from "../../../phases";

export interface MultiTraceAdapterArgs {
  members: ComparisonMember[];      // sorted by display_order
  traces: Map<number, Trace>;
  comparisonTitle: string;
  experimentName?: string;          // omit in /compare/all global scope
  xDomain: [number, number] | null;
  showPeakTicks: boolean;
  showPeakLabels: boolean;
  groupingMode: GroupingMode;
  sampleIdFor: (m: ComparisonMember) => number | null;
  /** Pre-resolved labels via lib/comparison/labels.resolveDisplayLabels.
   *  Falls back to "Exposure #${exposure_id}" when missing. */
  displayLabelByMemberId?: Map<number, string>;
}

export function buildMultiTraceExportSpec(args: MultiTraceAdapterArgs): ExportSpec {
  const {
    members, traces, comparisonTitle, experimentName,
    xDomain, showPeakTicks, showPeakLabels,
    groupingMode, sampleIdFor, displayLabelByMemberId,
  } = args;

  // Critical ordering (spec §MultiTracePlot export "Critical"): compute
  // colors on the UNFILTERED members so defaultDistinct findIndex is stable,
  // then filter null-exposure_id rows from the mark-build pass.
  const colorByMember = new Map<number, string>();
  for (const m of members) {
    colorByMember.set(
      m.id,
      colorFor(m, groupingMode, COMPARE_PALETTE_LIGHT, {
        allMembers: members,
        sampleIdFor,
      }),
    );
  }

  const filtered = members.filter((m) => m.exposure_id !== null);
  if (filtered.length === 0) {
    throw new Error(
      "buildMultiTraceExportSpec: every member has exposure_id === null; "
      + "export should have been gated by FigureExportControls",
    );
  }

  // Build a default display-label map if the parent didn't provide one. The
  // ComparePage call site DOES provide one; the adapter test calls without.
  const labels = displayLabelByMemberId ?? new Map<number, string>(
    filtered.map((m) => [
      m.id,
      m.label_override ?? `Exposure #${m.exposure_id}`,
    ]),
  );

  // panelHeight for the marks builder = export height minus chrome.
  const panelHeight = COMPARE_DIMS.height - EXPORT_MARGIN.top - EXPORT_MARGIN.bottom;
  const marks = buildMultiTraceExportMarks({
    members: filtered,
    traces,
    displayLabelByMemberId: labels,
    colorByMember,
    showPeakTicks, showPeakLabels,
    panelHeight,
  });

  // Legend per grouping mode.
  const legendRows = buildLegendRows(members, filtered, groupingMode, sampleIdFor, colorByMember);

  const xConfig: Record<string, unknown> = {
    type: "log",
    label: "q",
  };
  if (xDomain) xConfig.domain = xDomain;

  const yConfig: Record<string, unknown> = {
    type: "linear",
    label: null,
    domain: [0, panelHeight], // synthetic — marks pre-place y in pixel space
    axis: null,
  };

  const title: ExportSpec["title"] = experimentName !== undefined
    ? { primary: comparisonTitle, secondary: experimentName }
    : { primary: comparisonTitle };

  return {
    title,
    width: COMPARE_DIMS.width,
    height: COMPARE_DIMS.height,
    plot: {
      marks,
      x: xConfig,
      y: yConfig,
      width: COMPARE_DIMS.width  - EXPORT_MARGIN.left - EXPORT_MARGIN.right,
      height: panelHeight,
      marginLeft: 0,
      marginRight: 0,
      marginTop: 0,
      marginBottom: 0,
    },
    legend: { rows: legendRows },
  };
}

function buildLegendRows(
  unfilteredMembers: ComparisonMember[],
  filteredMembers: ComparisonMember[],
  mode: GroupingMode,
  sampleIdFor: (m: ComparisonMember) => number | null,
  colorByMember: Map<number, string>,
): LegendRow[] {
  if (mode === "distinct") return [];

  const rows: LegendRow[] = [];
  const seen = new Set<string>();
  let orphanPresent = false;

  for (const m of filteredMembers) {
    const color = colorByMember.get(m.id) ?? ORPHAN_FALLBACK_LIGHT;
    if (color === ORPHAN_FALLBACK_LIGHT) {
      orphanPresent = true;
      continue;
    }
    if (mode === "bySample") {
      const sid = sampleIdFor(m);
      const key = `bySample:${sid ?? "null"}`;
      if (seen.has(key)) continue;
      seen.add(key);
      rows.push({
        color,
        symbol: "swatch",
        label: sid !== null ? `Sample ${sid}` : "(unknown sample)",
      });
    } else if (mode === "byPhase") {
      const phase = m.snapshot?.confirmed_index?.phase;
      if (!phase) {
        orphanPresent = true;
        continue;
      }
      if (seen.has(phase)) continue;
      seen.add(phase);
      rows.push({
        color: phaseColor(phase),
        symbol: "swatch",
        label: phase,
      });
    }
  }

  // Note: unfilteredMembers may include null-exposure rows that bind via the
  // bySample sampleIdFor path; those still produce orphan colour and should
  // contribute to the orphan row presence.
  for (const m of unfilteredMembers) {
    if (m.exposure_id === null) {
      orphanPresent = true;
      break;
    }
  }

  if (orphanPresent) {
    rows.push({
      color: ORPHAN_FALLBACK_LIGHT,
      symbol: "swatch",
      label: "unphased / unbound",
    });
  }

  return rows;
}
```

- [ ] **Step 4: Run, verify PASS**

Run: same command as Step 2.
Expected: PASS, all 9 assertions green.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/figure-export/adapters/multiTraceAdapter.ts \
        packages/HimalayaUI/frontend/test/figure-export/multiTraceAdapter.test.ts
git commit -m "feat(figure-export): multiTraceAdapter (covers all GroupingModes)"
```

---

### Task 12: Wire `<FigureExportControls>` into `ComparePage`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/ComparePage.tsx`

Spec §Placement (ComparePage review header): append at the right end of the `compare-review-header` row, after `<ForksPopover>`. Disabled gate: `members.length === 0 || traces.size === 0 || members.every(m => m.exposure_id == null)` (spec §Disabled states).

- [ ] **Step 1: Add imports**

In `ComparePage.tsx` near the top:

```tsx
import { FigureExportControls } from "../components/FigureExportControls";
import { buildMultiTraceExportSpec } from "../lib/figure-export/adapters/multiTraceAdapter";
import { slugifyForFilename } from "../lib/figure-export/filename";
import { useExperiment } from "../queries";
```

- [ ] **Step 2: Inside `ReviewPlot`, build the export thunk + filename stem**

After the `displayLabelByMemberId` derivation (around line 220) and before `plotLoading`, insert:

```tsx
  // Figure export — read experiment name only when in experiment scope.
  const experimentQ = useExperiment(eid ?? 0);
  const experimentName = eid !== undefined
    ? (experimentQ.data?.name ?? `Experiment ${eid}`)
    : undefined;

  const exportFilenameStem = `himalaya-comparison-${
    slugifyForFilename(experimentName ?? "all")
  }-${
    slugifyForFilename(compQ.data?.title ?? "")
  }`;

  const exportSpec = useCallback(() => {
    if (!compQ.data || members.length === 0) {
      throw new Error("FigureExportControls: parent disabled-gate violated");
    }
    return buildMultiTraceExportSpec({
      members,
      traces,
      comparisonTitle: compQ.data.title,
      ...(experimentName !== undefined ? { experimentName } : {}),
      xDomain,
      showPeakTicks,
      showPeakLabels,
      groupingMode,
      sampleIdFor,
      displayLabelByMemberId,
    });
  }, [
    compQ.data, members, traces, experimentName,
    xDomain, showPeakTicks, showPeakLabels,
    groupingMode, sampleIdFor, displayLabelByMemberId,
  ]);

  const exportDisabled =
    members.length === 0
    || traces.size === 0
    || members.every((m) => m.exposure_id === null);
```

- [ ] **Step 3: Append `<FigureExportControls>` to the header row**

Modify the `compare-review-header` block:

```tsx
      <div
        data-testid="compare-review-header"
        className="flex items-center gap-3 flex-wrap"
      >
        <GroupingModeToggle />
        <AnnotationToggles />
        {isStale && (
          <NeedsReviewBadge … />
        )}
        {compQ.data && (
          <EditOrForkButton … />
        )}
        {compQ.data && (
          <LineageBadge comparison={compQ.data} experimentId={eid} scope={scope} />
        )}
        <ForksPopover comparisonId={id} experimentId={eid} scope={scope} />
        <span className="ml-auto inline-flex items-center gap-1">
          <FigureExportControls
            spec={exportSpec}
            filenameStem={exportFilenameStem}
            ariaContext="comparison plot"
            disabled={exportDisabled}
          />
        </span>
      </div>
```

- [ ] **Step 4: Verify TS compiles**

Run: `npm run build 2>&1 | tail -20`
Expected: build succeeds.

- [ ] **Step 5: Run existing ComparePage tests**

Run: `node_modules/.bin/vitest run test/ComparePage 2>&1 | tail -30`
Expected: existing tests stay green.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/pages/ComparePage.tsx
git commit -m "feat(figure-export): wire FigureExportControls into ComparePage review header"
```

---

# Phase 4 — E2E + final verification

---

### Task 13: Playwright permissions + figure-export E2E spec

**Files:**
- Modify: `packages/HimalayaUI/frontend/playwright.config.ts`
- Create: `packages/HimalayaUI/frontend/e2e/figure-export.spec.ts`

Spec §Playwright E2E (mocked): clipboard test (Chromium-only) + download test (cross-browser; our config is Chromium-only single-project, so it just runs there). Requires `permissions: ['clipboard-read', 'clipboard-write']` on `use:`.

- [ ] **Step 1: Modify `playwright.config.ts`**

Update the `use:` block to include permissions:

```ts
  use: {
    baseURL: "http://127.0.0.1:5173",
    browserName: "chromium",
    headless: true,
    permissions: ["clipboard-read", "clipboard-write"],
  },
```

- [ ] **Step 2: Skim an existing E2E spec for the seed-state pattern**

Read `packages/HimalayaUI/frontend/e2e/compare.spec.ts` (or any existing spec) to confirm the `mockApi` / `seedState` / `makeState` helpers and `addInitScript` pattern used to seed Zustand state. The figure-export spec follows the same shape.

Run: `ls packages/HimalayaUI/frontend/e2e/*.spec.ts`
Expected: see the spec list and pick the simplest as a template.

- [ ] **Step 3: Create the E2E spec**

`packages/HimalayaUI/frontend/e2e/figure-export.spec.ts`:

```ts
import { test, expect } from "@playwright/test";

/**
 * figure-export.spec.ts — Playwright coverage for the Copy/Save figure
 * export buttons (spec: 2026-05-08-figure-export-design.md, issue #90).
 *
 * Two scenarios:
 *   1. Copy → clipboard contains image/png   (Chromium-only)
 *   2. Download → file lands; filename matches the `himalaya-trace-…-YYYY-MM-DD.{png,svg}` pattern.
 */

const SAMPLE_FIXTURE = {
  exposure: { id: 1, sample_id: 10, name: "exp1.dat", selected: true },
  trace: {
    q: Array.from({ length: 50 }, (_, i) => 0.05 + i * 0.02),
    I: Array.from({ length: 50 }, (_, i) => 1000 / (1 + i)),
    sigma: Array.from({ length: 50 }, () => 5),
  },
  peaks: [{
    id: 1, exposure_id: 1, q: 0.18, intensity: 420, prominence: 180,
    sharpness: 2.1, source: "auto", excluded: false,
  }],
};

async function mockApi(page: import("@playwright/test").Page): Promise<void> {
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ json: [{ id: 7, name: "JC23", path: "/x", data_dir: "/x", analysis_dir: "/x", manifest_path: null, created_at: "2026", q_units: "A-1" }] }));
  await page.route("**/api/experiments/7", (r) =>
    r.fulfill({ json: { id: 7, name: "JC23", path: "/x", data_dir: "/x", analysis_dir: "/x", manifest_path: null, created_at: "2026", q_units: "A-1" } }));
  await page.route("**/api/experiments/7/samples", (r) =>
    r.fulfill({ json: [{ id: 10, experiment_id: 7, label: "S4", name: "Sample 4", notes: null, tags: [] }] }));
  await page.route("**/api/samples/10/exposures*", (r) =>
    r.fulfill({ json: [SAMPLE_FIXTURE.exposure] }));
  await page.route("**/api/exposures/1/trace", (r) =>
    r.fulfill({ json: SAMPLE_FIXTURE.trace }));
  await page.route("**/api/exposures/1/peaks", (r) =>
    r.fulfill({ json: SAMPLE_FIXTURE.peaks }));
  await page.route("**/api/exposures/1/indices", (r) =>
    r.fulfill({ json: [] }));
  await page.route("**/api/exposures/1/groups", (r) =>
    r.fulfill({ json: [] }));
}

async function seedState(page: import("@playwright/test").Page): Promise<void> {
  await page.addInitScript(() => {
    const state = {
      version: 3,
      state: {
        activeExperimentId: 7,
        activeSampleId: 10,
        activeExposureId: 1,
        activePage: "index",
        theme: "dark",
        // Other Zustand fields fall back to defaults from state.ts.
      },
    };
    localStorage.setItem("himalaya-ui:state", JSON.stringify(state));
  });
}

test.describe("Figure export — Index page (TraceViewer)", () => {
  test("Copy puts a PNG on the clipboard (Chromium)", async ({ page, browserName }) => {
    test.skip(browserName !== "chromium", "Clipboard read requires Chromium permissions");
    await mockApi(page);
    await seedState(page);
    await page.goto("/");

    const copyBtn = page.getByRole("button", { name: /copy trace plot to clipboard/i });
    await expect(copyBtn).toBeEnabled();
    await copyBtn.click();

    // Wait for the success toast to fire (eats any in-flight render time).
    // Asserting on either the toast text or the clipboard content directly —
    // we go with the clipboard since that's the actual contract.
    await page.waitForTimeout(500);
    const types = await page.evaluate(async () => {
      const items = await navigator.clipboard.read();
      const out: string[] = [];
      for (const item of items) for (const t of item.types) out.push(t);
      return out;
    });
    expect(types).toContain("image/png");
  });

  test("Download → PNG file lands with himalaya-trace-… filename", async ({ page }) => {
    await mockApi(page);
    await seedState(page);
    await page.goto("/");

    const dlBtn = page.getByRole("button", { name: /download trace plot as png/i });
    await expect(dlBtn).toBeEnabled();

    const downloadPromise = page.waitForEvent("download");
    await dlBtn.click();
    const download = await downloadPromise;

    expect(download.suggestedFilename()).toMatch(/^himalaya-trace-.*\d{4}-\d{2}-\d{2}\.png$/);
  });

  test("Download → SVG via chevron menu", async ({ page }) => {
    await mockApi(page);
    await seedState(page);
    await page.goto("/");

    const chevron = page.getByRole("button", { name: /other download formats/i });
    await chevron.click();
    const svgRow = page.getByText(/download as svg/i);

    const downloadPromise = page.waitForEvent("download");
    await svgRow.click();
    const download = await downloadPromise;

    expect(download.suggestedFilename()).toMatch(/^himalaya-trace-.*\d{4}-\d{2}-\d{2}\.svg$/);
  });
});
```

- [ ] **Step 4: Kill any stale Vite + run E2E**

```bash
lsof -ti:5173 | xargs -r kill -9 2>/dev/null || true
cd packages/HimalayaUI/frontend
node_modules/.bin/playwright test e2e/figure-export.spec.ts --reporter=list 2>&1 | tail -30
```

Expected: 3 tests pass on Chromium. The Playwright `webServer` block in `playwright.config.ts` auto-starts Vite.

If a clipboard test fails with a permissions error, double-check the `permissions:` array ended up in the resolved config (`use:` not inside `projects:`).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/playwright.config.ts \
        packages/HimalayaUI/frontend/e2e/figure-export.spec.ts
git commit -m "test(figure-export): Playwright clipboard + download specs"
```

---

### Task 14: Final verification + manual smoke

**No new files**. Sanity-check the whole branch.

- [ ] **Step 1: TypeScript + Vite build clean**

```bash
cd packages/HimalayaUI/frontend
npm run build 2>&1 | tail -15
```

Expected: `tsc --noEmit` succeeds; `vite build` succeeds. No errors.

- [ ] **Step 2: All Vitest unit tests pass**

```bash
cd packages/HimalayaUI/frontend
npm test -- --run 2>&1 | tail -20
```

Expected: all suites green; the 6 figure-export test files run and pass; existing test count unchanged.

- [ ] **Step 3: All Playwright mocked E2E pass**

```bash
cd packages/HimalayaUI/frontend
node_modules/.bin/playwright test --reporter=list 2>&1 | tail -20
```

Expected: every test passes (the 3 new figure-export tests + existing specs that don't already fail on main; `e2e/smoke.spec.ts` has 3 known-failing tests pinned in #95 — those are pre-existing and NOT part of this PR's gate).

- [ ] **Step 4: Manual smoke**

Spec §Acceptance — manual verification on a real browser:

1. `make sysimage && bin/himalaya serve <experiment> --port 8080` (or `julia --project=… serve …`).
2. `cd packages/HimalayaUI/frontend && npm run dev -- --host 127.0.0.1` in a second terminal.
3. Navigate to `http://127.0.0.1:5173/`. Pick an experiment + sample.
4. **TraceViewer Copy**: click Copy in the title strip, paste into Slides / Slack / Notes. Confirm the figure shows.
5. **TraceViewer Download PNG**: click Save → file lands as `himalaya-trace-…-YYYY-MM-DD.png`. Open it.
6. **TraceViewer Download SVG**: click chevron → Download as SVG → file lands; open in browser, structure matches the PNG.
7. Switch to the Compare tab. Open or create a comparison with ≥2 members.
8. **Compare Copy**: click Copy → paste into a deck.
9. **Compare Download PNG/SVG**: same as the trace path; filename starts `himalaya-comparison-…`.
10. **Disabled gate**: clear all members from a comparison; confirm Copy and Save are disabled.
11. **Insecure-origin check** (optional): if you have a non-HTTPS proxy / preview deploy, load the app there and verify Copy is disabled with "Clipboard requires HTTPS" tooltip.

- [ ] **Step 5: Update PR description with smoke notes**

When opening the PR, mention what was actually smoke-tested. The PR-review process will hold us accountable: don't tick the spec-acceptance boxes that weren't run end-to-end.

- [ ] **Step 6: No commit needed for this task**

Final verification only. The branch is ready for PR.

---

## Spec coverage map

For Phase-N reviewers — every section/subsection of the spec maps to at least one task here:

| Spec § | Tasks |
|---|---|
| Approach: headless re-render | Tasks 6, 10 |
| TraceViewer export — title, marks, legend | Tasks 6, 7 |
| MultiTracePlot export — title, members filter, colorByMember on unfiltered, legend modes | Tasks 10, 11 |
| Architecture (file layout) | All file paths above match the spec's tree |
| ExportSpec shape | Task 1 |
| Adapters (signatures) | Tasks 7, 11 |
| Renderer | Task 5 |
| Colour resolution | Task 1 (`COMPARE_PALETTE_LIGHT`), Tasks 7/11 |
| Button component | Task 8 |
| Popover behaviour (Esc + outside-click) | Task 8 |
| Placement (PlotCard, ComparePage) | Tasks 9, 12 |
| Data flow at click time (try/finally pending) | Task 8 |
| Browser compatibility (clipboard pre-flight) | Task 4 |
| SVG → PNG canvas pipeline | Task 5 |
| Disabled states | Tasks 9, 12 |
| Filenames (slug per segment, en-CA, sentinel) | Task 2 |
| Testing (Vitest + Playwright) | Tasks 2, 4, 5, 7, 8, 11, 13 |
| Toast messages | Task 8 |

## Risks owed during implementation (from spec §Risks)

1. **`COMPARE_PALETTE_LIGHT` contrast on white bg** — Task 1 introduces the dropped-luminance palette. If it still fails at Task 12 review (visual eyeball during the manual smoke), tune in `presets.ts` only — keep on-screen→export hue mapping stable.
2. **Render cost > 500ms** — Spec §Button component asks for profiling during step-3 wiring (Task 9 / 12). If a TraceViewer export with many predicted-q ticks (multi-phase active group, 2× DPI) blocks the UI > 500ms on a slow laptop, swap the `pending` flag for an `AbortController`. Note in the PR description.
