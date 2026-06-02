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

test("LUT is non-inverting — brighter source pixel -> lighter output", async () => {
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
    expect(d[0]).toBeGreaterThanOrEqual(d[2]);       // warm at the dark end (R >= B)
  } finally {
    vi.restoreAllMocks();
  }
});

/**
 * U-3: drive the orient decision toward landscape (wide container, tall
 * viewport, square image) and assert a THUMB stays locked to portrait while a
 * FULL frame still rotates. JSDOM has no layout, so clientWidth/clientHeight
 * are stubbed via getters; window.innerWidth is forced past ROTATE_MIN_VIEWPORT.
 */
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
  global.createImageBitmap = vi.fn().mockResolvedValue({ width: 8, height: 8, close: vi.fn() } as unknown as ImageBitmap);
  const restore = forceWideGeometry();
  try {
    render(<DetectorImage src="/x.png" size="thumb" />);
    const wrapper = await waitFor(() => {
      const el = screen.getByRole("img", { hidden: true }).closest("[data-orient]") as HTMLElement;
      expect(el).toBeTruthy();
      return el;
    });
    expect(wrapper).toHaveAttribute("data-orient", "portrait");
    const canvas = screen.getByRole("img", { hidden: true }) as HTMLCanvasElement;
    expect(canvas.style.transform).toBe("");
  } finally {
    restore();
  }
});

test("U-3 regression: a full frame still rotates under the same wide geometry", async () => {
  global.createImageBitmap = vi.fn().mockResolvedValue({ width: 8, height: 8, close: vi.fn() } as unknown as ImageBitmap);
  const restore = forceWideGeometry();
  try {
    render(<DetectorImage src="/x.png" size="full" />);
    await waitFor(() => {
      const canvas = screen.getByRole("img", { hidden: true }) as HTMLCanvasElement;
      expect(canvas.style.transform).toContain("rotate(90deg)");
    });
  } finally {
    restore();
  }
});
