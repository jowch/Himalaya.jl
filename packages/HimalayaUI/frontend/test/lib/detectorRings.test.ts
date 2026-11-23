import { describe, it, expect } from "vitest";
import { detectorGeometry, ringCircles } from "../../src/lib/detectorRings";
import type { Setup } from "../../src/api";

const setup = (over: Partial<Setup> = {}): Setup => ({
  beam_center_x: 421.409, beam_center_y: 836.946, pixel_size_um: 172, energy_kev: 9.000027604502573, flight_path_m: 1.8095, q_units: "A-1",
  ...over,
} as Setup);

describe("detectorGeometry", () => {
  it("needs the beam centre, pixel size, energy and flight path", () => {
    expect(detectorGeometry(setup())).not.toBeNull();
    for (const k of ["beam_center_x", "beam_center_y", "pixel_size_um", "energy_kev", "flight_path_m"] as const) {
      expect(detectorGeometry(setup({ [k]: null }))).toBeNull();
      expect(detectorGeometry(setup({ [k]: 0 }))).toBeNull();
    }
    expect(detectorGeometry(undefined)).toBeNull();
  });
});

describe("ringCircles", () => {
  const size = { width: 981, height: 1043 };

  it("puts the beam centre in image pixels, measured from the bottom-left with y up", () => {
    const c = ringCircles([0.1], detectorGeometry(setup())!, size)!;
    expect(c.cx).toBeCloseTo(421.909, 3);
    expect(c.cy).toBeCloseTo(1043 - 836.946 - 0.5, 3);
  });

  it("matches the radius the demo generator draws a ring at, q = 4π sin(θ)/λ with r = L tan 2θ", () => {
    const lambda = 12398.42 / 9000.027604502573;
    const radius = (q: number) => (1.8095 * Math.tan(2 * Math.asin((q * lambda) / (4 * Math.PI)))) / 172e-6;
    const c = ringCircles([0.05, 0.1, 0.2], detectorGeometry(setup())!, size)!;
    c.radii.forEach((r, i) => expect(r).toBeCloseTo(radius([0.05, 0.1, 0.2][i]!), 1));
    expect(c.radii[1]).toBeCloseTo(230.7, 1);
  });

  it("reads q in nm⁻¹ as ten times q in Å⁻¹", () => {
    const a = ringCircles([0.1], detectorGeometry(setup())!, size)!;
    const nm = ringCircles([1], detectorGeometry(setup({ q_units: "nm-1" }))!, size)!;
    expect(nm.radii[0]).toBeCloseTo(a.radii[0]!, 6);
  });

  it("skips orders past backscatter and returns nothing without a size", () => {
    const c = ringCircles([0.1, 100], detectorGeometry(setup())!, size)!;
    expect(c.radii).toHaveLength(1);
    expect(ringCircles([0.1], detectorGeometry(setup())!, undefined)).toBeNull();
  });
});
