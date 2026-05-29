import { render, screen, waitFor } from "@testing-library/react";
import { vi, beforeEach } from "vitest";
import { DetectorImage } from "../src/components/DetectorImage";

// Minimal 1×1 white PNG (base64)
const TINY_PNG =
  "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwADhQGAWjR9awAAAABJRU5ErkJggg==";

beforeEach(() => {
  global.fetch = vi.fn().mockResolvedValue({
    ok: true,
    blob: () =>
      Promise.resolve(
        new Blob(
          [Uint8Array.from(atob(TINY_PNG), (c) => c.charCodeAt(0))],
          { type: "image/png" },
        ),
      ),
  } as Response);

  global.createImageBitmap = vi.fn().mockResolvedValue({
    width: 1,
    height: 1,
    close: vi.fn(),
  } as unknown as ImageBitmap);

  // OffscreenCanvas is not available in JSDOM
  const mockOffscreen = {
    getContext: () => ({
      drawImage: vi.fn(),
      getImageData: () => ({ data: new Uint8ClampedArray(4) }),
    }),
  };
  // @ts-expect-error JSDOM stub
  global.OffscreenCanvas = vi.fn().mockImplementation(() => mockOffscreen);
});

test("canvas has correct dimensions after bitmap.close() neuters width/height", async () => {
  // Real ImageBitmap.close() sets width and height to 0 (spec: "neutered").
  // This test ensures we capture dimensions BEFORE calling close().
  let bitmapClosed = false;
  const mockBitmap = Object.defineProperties({} as ImageBitmap, {
    width:  { get: () => (bitmapClosed ? 0 : 8) },
    height: { get: () => (bitmapClosed ? 0 : 6) },
    close:  { value: vi.fn(() => { bitmapClosed = true; }) },
  });
  global.createImageBitmap = vi.fn().mockResolvedValue(mockBitmap);

  render(
    <DetectorImage exposureId={1} imagePath="/tmp/foo.tiff"
      imageVersion="v1-1700000000" size="full" />,
  );

  await waitFor(() => {
    const canvas = screen.getByRole("img", { hidden: true }) as HTMLCanvasElement;
    // If the bug were present (reading bitmap.width after close), canvas
    // would be 0×0 and these assertions would fail.
    expect(canvas.width).toBe(8);
    expect(canvas.height).toBe(6);
  });
});

test("renders a canvas element when imagePath is provided", async () => {
  render(
    <DetectorImage exposureId={1} imagePath="/tmp/foo.tiff"
      imageVersion="v1-1700000000" size="full" />,
  );
  await waitFor(() =>
    expect(screen.getByRole("img", { hidden: true })).toBeInTheDocument(),
  );
});

test("shows placeholder when imagePath is null", () => {
  render(<DetectorImage exposureId={1} imagePath={null}
    imageVersion="" size="full" />);
  expect(
    screen.getByTestId("detector-image-placeholder"),
  ).toBeInTheDocument();
});

test("appends ?v=<imageVersion> and does not request `cache: no-store`", async () => {
  // Versioned URL allows browser caching while staying correct: a new mtime
  // or processing-version bump produces a new URL → fresh fetch.
  const fetchSpy = vi.fn().mockResolvedValue({
    ok: true,
    blob: () => Promise.resolve(new Blob([new Uint8Array(0)], { type: "image/png" })),
  } as Response);
  global.fetch = fetchSpy;

  render(
    <DetectorImage exposureId={42} imagePath="/tmp/x.tiff"
      imageVersion="v1-1700000099" size="full" />,
  );

  await waitFor(() => expect(fetchSpy).toHaveBeenCalled());
  const [url, init] = fetchSpy.mock.calls[0];
  expect(url).toBe("/api/exposures/42/image?v=v1-1700000099");
  // Default fetch (no `cache: "no-store"`) so the browser cache can serve
  // stable URLs across theme toggles and selection changes.
  expect(init).toBeUndefined();
});

test("thumb URL preserves `thumb=1` alongside the version param", async () => {
  const fetchSpy = vi.fn().mockResolvedValue({
    ok: true,
    blob: () => Promise.resolve(new Blob([new Uint8Array(0)], { type: "image/png" })),
  } as Response);
  global.fetch = fetchSpy;

  render(
    <DetectorImage exposureId={7} imagePath="/tmp/x.tiff"
      imageVersion="v1-42" size="thumb" />,
  );

  await waitFor(() => expect(fetchSpy).toHaveBeenCalled());
  expect(fetchSpy.mock.calls[0][0]).toBe("/api/exposures/7/image?thumb=1&v=v1-42");
});

/**
 * U-3 (#256): drive the orient decision toward landscape (wide container, tall
 * viewport, square image) and assert a THUMB stays locked to portrait while a
 * FULL frame still rotates. JSDOM has no layout, so clientWidth/clientHeight
 * are stubbed via getters; window.innerWidth is forced past ROTATE_MIN_VIEWPORT.
 */
function forceWideGeometry(): () => void {
  const protoW = Object.getOwnPropertyDescriptor(
    HTMLElement.prototype, "clientWidth",
  );
  const protoH = Object.getOwnPropertyDescriptor(
    HTMLElement.prototype, "clientHeight",
  );
  Object.defineProperty(HTMLElement.prototype, "clientWidth", {
    configurable: true, get() { return 400; },
  });
  Object.defineProperty(HTMLElement.prototype, "clientHeight", {
    configurable: true, get() { return 40; },
  });
  const origInner = window.innerWidth;
  Object.defineProperty(window, "innerWidth", {
    configurable: true, value: 1600,
  });
  return () => {
    if (protoW) Object.defineProperty(HTMLElement.prototype, "clientWidth", protoW);
    if (protoH) Object.defineProperty(HTMLElement.prototype, "clientHeight", protoH);
    Object.defineProperty(window, "innerWidth", {
      configurable: true, value: origInner,
    });
  };
}

test("U-3: a thumb stays portrait even when geometry would rotate a full frame", async () => {
  // 8×8 (square) image, container 400×40 (10:1) → decideOrient would pick
  // landscape for a full frame. A thumb must IGNORE that and lock portrait.
  global.createImageBitmap = vi.fn().mockResolvedValue({
    width: 8, height: 8, close: vi.fn(),
  } as unknown as ImageBitmap);
  const restore = forceWideGeometry();
  try {
    render(
      <DetectorImage exposureId={9} imagePath="/tmp/x.tiff"
        imageVersion="v1-9" size="thumb" />,
    );
    const wrapper = await waitFor(() => {
      const el = screen
        .getByRole("img", { hidden: true })
        .closest("[data-orient]") as HTMLElement;
      expect(el).toBeTruthy();
      return el;
    });
    // The gate holds: a thumb never flips to landscape.
    expect(wrapper).toHaveAttribute("data-orient", "portrait");
    const canvas = screen.getByRole("img", { hidden: true }) as HTMLCanvasElement;
    expect(canvas.style.transform).toBe("");
  } finally {
    restore();
  }
});

test("U-3 regression: a full frame still rotates under the same wide geometry", async () => {
  // Same geometry as the thumb test, but size="full" → the auto-rotate path is
  // untouched, so the canvas carries a rotate transform. Proves the gate is
  // size-scoped, not a blanket disable.
  global.createImageBitmap = vi.fn().mockResolvedValue({
    width: 8, height: 8, close: vi.fn(),
  } as unknown as ImageBitmap);
  const restore = forceWideGeometry();
  try {
    render(
      <DetectorImage exposureId={9} imagePath="/tmp/x.tiff"
        imageVersion="v1-9" size="full" />,
    );
    await waitFor(() => {
      const canvas = screen.getByRole("img", { hidden: true }) as HTMLCanvasElement;
      expect(canvas.style.transform).toContain("rotate(90deg)");
    });
  } finally {
    restore();
  }
});

/**
 * B.1 / U-1 (#255): the detector LUT must warm the image to the detector-window
 * tokens AND stay NON-INVERTING — the brightest input pixel maps to the LIGHTER
 * (signal) output endpoint, the darkest to the near-black window backing. The
 * sign of that ramp is the actual perceptual fix (saxs test-pin) and the thing
 * most likely to silently regress.
 */
test("B.1: LUT is non-inverting — brighter intensity → lighter output", async () => {
  // Two source pixels: a dark one (intensity 0) and a bright one (intensity 255).
  // The component reads the R channel as t and rebuilds RGB in place; we capture
  // the buffer handed to putImageData.
  const srcData = new Uint8ClampedArray([
    0, 0, 0, 255, // intensity 0  → window backing (darkest)
    255, 255, 255, 255, // intensity 255 → signal (lightest)
  ]);
  const mockOffscreen = {
    getContext: () => ({
      drawImage: vi.fn(),
      getImageData: () => ({ data: srcData }),
    }),
  };
  // @ts-expect-error JSDOM stub
  global.OffscreenCanvas = vi.fn().mockImplementation(() => mockOffscreen);

  // Stub getComputedStyle so getCssColor resolves the two warm endpoints to
  // known sRGB triples without depending on JSDOM oklch support. frame-edge is
  // near-black; frame-signal is warm off-white.
  const realGCS = window.getComputedStyle;
  vi.spyOn(window, "getComputedStyle").mockImplementation(((el: Element) => {
    if (el === document.documentElement) {
      return {
        getPropertyValue: (name: string) =>
          name === "--color-frame-edge"
            ? "rgb(20, 18, 14)"
            : name === "--color-frame-signal"
              ? "rgb(238, 236, 228)"
              : "",
      } as CSSStyleDeclaration;
    }
    return realGCS(el);
  }) as typeof window.getComputedStyle);

  global.createImageBitmap = vi.fn().mockResolvedValue({
    width: 2, height: 1, close: vi.fn(),
  } as unknown as ImageBitmap);

  let captured: Uint8ClampedArray | null = null;
  // getCssColor() creates its OWN 1x1 canvas, sets fillStyle to the resolved CSS
  // string, fills, then reads back the pixel to parse it to RGB. JSDOM's 2d
  // context lacks these, so the stub doubles as a tiny rgb() parser: getImageData
  // returns the channels from the last `rgb(r, g, b)` fillStyle. The on-screen
  // canvas's putImageData is captured for the LUT assertion.
  vi.spyOn(HTMLCanvasElement.prototype, "getContext").mockImplementation(function () {
    let lastFill = "rgb(0, 0, 0)";
    return {
      set fillStyle(v: string) { lastFill = v; },
      get fillStyle() { return lastFill; },
      fillRect: () => {},
      drawImage: () => {},
      getImageData: () => {
        const m = /rgba?\((\d+),\s*(\d+),\s*(\d+)/.exec(lastFill);
        const [r, g, b] = m ? [+m[1], +m[2], +m[3]] : [0, 0, 0];
        return { data: new Uint8ClampedArray([r, g, b, 255]) };
      },
      putImageData: (imgData: ImageData) => {
        captured = imgData.data;
      },
    } as unknown as CanvasRenderingContext2D;
  } as typeof HTMLCanvasElement.prototype.getContext);

  try {
    render(
      <DetectorImage exposureId={1} imagePath="/tmp/x.tiff"
        imageVersion="v3-1" size="full" />,
    );
    await waitFor(() => expect(captured).not.toBeNull());
    const d = captured as unknown as Uint8ClampedArray;
    const darkSum = d[0] + d[1] + d[2];
    const brightSum = d[4] + d[5] + d[6];
    // Non-inverting: the bright input is lighter than the dark input.
    expect(brightSum).toBeGreaterThan(darkSum);
    // Warm endpoints: t=0 lands on frame-edge (near-black), t=255 on frame-signal.
    expect(d[0]).toBe(20); expect(d[1]).toBe(18); expect(d[2]).toBe(14);
    expect(d[4]).toBe(238); expect(d[5]).toBe(236); expect(d[6]).toBe(228);
  } finally {
    vi.restoreAllMocks();
  }
});

/**
 * B.4 / U-2 / R3-S06 (#255): the missing-image placeholder renders as a
 * `frame-edge` window with a `frame-tag` mono caption — not light `text-fg-muted`
 * text on paper. data-* + class-list assertions (no class-string matching of the
 * positive utilities; we assert the dark-era survivor is GONE).
 */
test("B.4: missing-image placeholder declares the frame-edge window variant", () => {
  render(<DetectorImage exposureId={1} imagePath={null}
    imageVersion="" size="full" />);
  const ph = screen.getByTestId("detector-image-placeholder");
  // `data-variant="frame-window"` is the behavioral contract that the empty
  // state renders as a frame-edge window (matching the live detector), not
  // light text on paper. Per test/AGENTS.md we assert the data-attribute, not
  // the Tailwind class strings that implement the window/caption styling
  // (R3-S06 — the legacy text-fg-muted treatment — is retired by that variant).
  expect(ph).toHaveAttribute("data-variant", "frame-window");
});
