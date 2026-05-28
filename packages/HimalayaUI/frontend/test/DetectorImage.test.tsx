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
