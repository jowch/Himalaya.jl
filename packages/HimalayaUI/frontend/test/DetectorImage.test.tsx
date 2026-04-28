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
    <DetectorImage exposureId={1} imagePath="/tmp/foo.tiff" size="full" />,
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
    <DetectorImage exposureId={1} imagePath="/tmp/foo.tiff" size="full" />,
  );
  await waitFor(() =>
    expect(screen.getByRole("img", { hidden: true })).toBeInTheDocument(),
  );
});

test("shows placeholder when imagePath is null", () => {
  render(<DetectorImage exposureId={1} imagePath={null} size="full" />);
  expect(
    screen.getByTestId("detector-image-placeholder"),
  ).toBeInTheDocument();
});
