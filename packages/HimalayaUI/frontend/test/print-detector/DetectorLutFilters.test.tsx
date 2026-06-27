import { render } from "@testing-library/react";
import { test, expect } from "vitest";
import { DetectorLutFilters } from "../../src/print/detector/DetectorLutFilters";
import { detectorLutTableValues } from "../../src/print/detector/detectorLut";

// Guards the provider wiring DetectorImage references via `filter: url(#…)`. The
// moved "LUT non-inverting" test covers the table MATH; this covers that the
// filters actually exist, carry all three channel funcs, and feed each channel
// its matching ramp (a swapped R/B or wrong #id would otherwise pass silently).
test("mounts both colormap filters with sRGB and all three channel ramps", () => {
  const { container } = render(<DetectorLutFilters />);
  for (const variant of ["neutral", "warm"] as const) {
    const filter = container.querySelector(`#detector-lut-${variant}`);
    expect(filter, `#detector-lut-${variant} present`).not.toBeNull();
    // sRGB is load-bearing (default linearRGB would gamma-shift the LUT).
    expect(filter!.getAttribute("color-interpolation-filters")).toBe("sRGB");

    const expected = detectorLutTableValues(variant);
    const funcs = filter!.querySelectorAll("feComponentTransfer > *");
    const byTag = (t: string) => filter!.querySelector(t);
    expect(funcs.length).toBe(3);
    expect(byTag("feFuncR")!.getAttribute("tableValues")).toBe(expected.r);
    expect(byTag("feFuncG")!.getAttribute("tableValues")).toBe(expected.g);
    expect(byTag("feFuncB")!.getAttribute("tableValues")).toBe(expected.b);
  }
});
