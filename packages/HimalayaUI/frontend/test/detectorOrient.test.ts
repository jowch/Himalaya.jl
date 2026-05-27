import { describe, it, expect } from "vitest";
import {
  decideOrient, ROTATE_THRESHOLD, ROTATE_MIN_VIEWPORT,
} from "../src/lib/detectorOrient";

describe("decideOrient", () => {
  it("forces portrait below the viewport gate", () => {
    const r = decideOrient({
      containerW: 800, containerH: 300, imageW: 100, imageH: 100,
      viewportW: ROTATE_MIN_VIEWPORT - 1,
    });
    expect(r.orient).toBe("portrait");
    expect(r.caps).toBeNull();
  });

  it("rotates to landscape when container is much wider than the image", () => {
    const r = decideOrient({
      containerW: 800, containerH: 300, imageW: 100, imageH: 100,
      viewportW: 1600,
    });
    expect(r.orient).toBe("landscape");
    // pre-rotation caps: maxW from container HEIGHT, maxH from container WIDTH
    expect(r.caps).toEqual({ maxW: 300, maxH: 800 });
  });

  it("stays portrait when aspect is within threshold", () => {
    const r = decideOrient({
      containerW: 110, containerH: 100, imageW: 100, imageH: 100,
      viewportW: 1600,
    });
    expect(110 / 100).toBeLessThan(ROTATE_THRESHOLD * (100 / 100));
    expect(r.orient).toBe("portrait");
  });

  it("returns portrait for degenerate zero sizes", () => {
    expect(decideOrient({
      containerW: 0, containerH: 0, imageW: 0, imageH: 0, viewportW: 1600,
    }).orient).toBe("portrait");
  });
});
