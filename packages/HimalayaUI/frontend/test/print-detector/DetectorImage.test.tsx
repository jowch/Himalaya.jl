import { render, screen, waitFor } from "@testing-library/react";
import { vi, beforeEach, test, expect } from "vitest";
import { DetectorImage } from "../../src/print/detector/DetectorImage";
import { detectorLutTableValues } from "../../src/print/detector/detectorLut";

const TINY_PNG =
  "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwADhQGAWjR9awAAAABJRU5ErkJggg==";

function pngBlob(): Blob {
  return new Blob([Uint8Array.from(atob(TINY_PNG), (c) => c.charCodeAt(0))], { type: "image/png" });
}

beforeEach(() => {
  global.fetch = vi.fn().mockResolvedValue({
    ok: true,
    headers: { get: (k: string) => (k === "X-Image-Width" ? "2048" : k === "X-Image-Height" ? "1024" : null) },
    blob: () => Promise.resolve(pngBlob()),
  } as unknown as Response);
  // JSDOM doesn't implement object URLs — stub them. The component creates one
  // per blob and revokes on swap/unmount; assert nothing throws + img gets a src.
  global.URL.createObjectURL = vi.fn(() => "blob:mock-url");
  global.URL.revokeObjectURL = vi.fn();
});

test("shows the frame-window placeholder when src is null", () => {
  render(<DetectorImage src={null} size="full" />);
  expect(screen.getByTestId("detector-image-placeholder")).toHaveAttribute("data-variant", "frame-window");
});

test("renders an <img> (role=img) fed by the fetched blob's object URL", async () => {
  render(<DetectorImage src="/fixtures/thumbs/37.png" size="full" />);
  const img = await waitFor(() => screen.getByRole("img", { hidden: true }) as HTMLImageElement);
  expect(img.tagName).toBe("IMG");
  await waitFor(() => expect(img.getAttribute("src")).toBe("blob:mock-url"));
});

test("fetches exactly the src it was given (no URL building inside)", async () => {
  const spy = vi.fn().mockResolvedValue({
    ok: true, blob: () => Promise.resolve(pngBlob()),
  } as Response);
  global.fetch = spy;
  render(<DetectorImage src="/api/exposures/42/image?thumb=1&v=v1-9" size="thumb" />);
  await waitFor(() => expect(spy).toHaveBeenCalled());
  expect(spy.mock.calls[0][0]).toBe("/api/exposures/42/image?thumb=1&v=v1-9");
  expect(spy.mock.calls[0][1]).toBeUndefined();
});

test("portrait image fills the frame (object-fit contain, scales up AND down)", async () => {
  render(<DetectorImage src="/x.png" size="full" />);
  const img = await waitFor(() => screen.getByRole("img", { hidden: true }) as HTMLImageElement);
  // Fill the frame, aspect-preserved — NOT maxWidth:100% which only scaled down
  // and left a sub-frame image at native pixel size.
  expect(img.style.objectFit).toBe("contain");
  expect(img.style.width).toBe("100%");
  expect(img.style.height).toBe("100%");
  expect(img.style.transform).toBe(""); // portrait: no rotate
});

test("recolors via the SVG colormap filter — neutral default, warm on request", async () => {
  const { rerender } = render(<DetectorImage src="/x.png" size="full" />);
  const img = await waitFor(() => screen.getByRole("img", { hidden: true }) as HTMLImageElement);
  expect(img.style.filter).toBe("url(#detector-lut-neutral)");
  // Switching variants restyles the SAME element — and must NOT refetch.
  const fetchCalls = (global.fetch as ReturnType<typeof vi.fn>).mock.calls.length;
  rerender(<DetectorImage src="/x.png" size="full" lutVariant="warm" />);
  expect(img.style.filter).toBe("url(#detector-lut-warm)");
  expect((global.fetch as ReturnType<typeof vi.fn>).mock.calls.length).toBe(fetchCalls);
});

test("LUT table values are non-inverting and warm at the dark end", () => {
  // The colormap correctness that the old canvas-loop test asserted now lives in
  // the SVG-filter table (one source: buildPrintDetectorLut). Verify the property
  // at the source rather than through a JSDOM canvas it can't actually paint.
  const { r, g, b } = detectorLutTableValues("neutral");
  const R = r.split(" ").map(Number), G = g.split(" ").map(Number), B = b.split(" ").map(Number);
  expect(R).toHaveLength(256);
  const lum = (i: number): number => R[i] + G[i] + B[i];
  expect(lum(255)).toBeGreaterThan(lum(0));   // brighter source → lighter output
  expect(R[0]).toBeGreaterThanOrEqual(B[0]);  // warm (R >= B) at the dark end
});

/**
 * U-3: drive the orient decision toward landscape (wide container, tall
 * viewport, square image) and assert a THUMB stays locked to portrait while a
 * FULL frame still rotates. Image dims now come from the X-Image-Width/Height
 * headers (aspect is preserved by the server's 1536 cap), so the tests feed a
 * square detector via the headers. JSDOM has no layout, so clientWidth/Height
 * are stubbed via getters; window.innerWidth is forced past ROTATE_MIN_VIEWPORT.
 */
function squareHeaders(): void {
  global.fetch = vi.fn().mockResolvedValue({
    ok: true,
    headers: { get: (k: string) => (k === "X-Image-Width" || k === "X-Image-Height" ? "800" : null) },
    blob: () => Promise.resolve(pngBlob()),
  } as unknown as Response);
}

function forceWideGeometry(): () => void {
  const protoW = Object.getOwnPropertyDescriptor(HTMLElement.prototype, "clientWidth");
  const protoH = Object.getOwnPropertyDescriptor(HTMLElement.prototype, "clientHeight");
  Object.defineProperty(HTMLElement.prototype, "clientWidth", {
    configurable: true, get() { return 400; },
  });
  Object.defineProperty(HTMLElement.prototype, "clientHeight", {
    configurable: true, get() { return 40; },
  });
  const origInner = window.innerWidth;
  Object.defineProperty(window, "innerWidth", { configurable: true, value: 1600 });
  return () => {
    if (protoW) Object.defineProperty(HTMLElement.prototype, "clientWidth", protoW);
    if (protoH) Object.defineProperty(HTMLElement.prototype, "clientHeight", protoH);
    Object.defineProperty(window, "innerWidth", { configurable: true, value: origInner });
  };
}

test("U-3: a thumb stays portrait even when geometry would rotate a full frame", async () => {
  squareHeaders();
  const restore = forceWideGeometry();
  try {
    render(<DetectorImage src="/x.png" size="thumb" />);
    const wrapper = await waitFor(() => {
      const el = screen.getByRole("img", { hidden: true }).closest("[data-orient]") as HTMLElement;
      expect(el).toBeTruthy();
      return el;
    });
    expect(wrapper).toHaveAttribute("data-orient", "portrait");
    const img = screen.getByRole("img", { hidden: true }) as HTMLImageElement;
    expect(img.style.transform).toBe("");
  } finally {
    restore();
  }
});

test("U-3 regression: a full frame still rotates under the same wide geometry", async () => {
  squareHeaders();
  const restore = forceWideGeometry();
  try {
    render(<DetectorImage src="/x.png" size="full" />);
    await waitFor(() => {
      const img = screen.getByRole("img", { hidden: true }) as HTMLImageElement;
      expect(img.style.transform).toContain("rotate(90deg)");
    });
  } finally {
    restore();
  }
});

test("reports raw image size from X-Image-Width/Height headers", async () => {
  const onRawSize = vi.fn();
  render(<DetectorImage src="/x.png" size="full" onRawSize={onRawSize} />);
  await waitFor(() => expect(onRawSize).toHaveBeenCalledWith(2048, 1024));
});

test("a headerless response does not call onRawSize and does not throw", async () => {
  global.fetch = vi.fn().mockResolvedValue({
    ok: true,
    blob: () => Promise.resolve(pngBlob()),
  } as unknown as Response);
  const onRawSize = vi.fn();
  render(<DetectorImage src="/x.png" size="full" onRawSize={onRawSize} />);
  await waitFor(() => expect(screen.getByRole("img", { hidden: true })).toBeInTheDocument());
  expect(onRawSize).not.toHaveBeenCalled();
});

test("headers present but missing the keys (get→null) does not call onRawSize", async () => {
  // Exercises the `Number(null) === 0` branch of the dual guard (>0), distinct
  // from the no-headers-object branch above.
  global.fetch = vi.fn().mockResolvedValue({
    ok: true,
    headers: { get: () => null },
    blob: () => Promise.resolve(pngBlob()),
  } as unknown as Response);
  const onRawSize = vi.fn();
  render(<DetectorImage src="/x.png" size="full" onRawSize={onRawSize} />);
  await waitFor(() => expect(screen.getByRole("img", { hidden: true })).toBeInTheDocument());
  expect(onRawSize).not.toHaveBeenCalled();
});
