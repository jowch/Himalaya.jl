import { describe, expect, it } from "vitest";
import { exposureKindLabel, exposureTimeLabel } from "../../src/workspace/exposureLabels";
import { makeExposure } from "../helpers/exposure";

describe("exposureKindLabel", () => {
  it("names a single image, an average and a subtraction", () => {
    expect(exposureKindLabel(makeExposure())).toBe("Single");
    expect(exposureKindLabel(makeExposure({ kind: "averaged" }))).toBe("Averaged");
    expect(exposureKindLabel(makeExposure({ kind: "background_subtracted" }))).toBe("Subtracted");
  });

  it("says how many images the beamline kept when it rejected some", () => {
    expect(exposureKindLabel(makeExposure({ frames_total: 10, frames_rejected: 1 })))
      .toBe("Beamline average of 9 of 10 images");
  });

  it("says how many images the beamline averaged when it rejected none", () => {
    expect(exposureKindLabel(makeExposure({ frames_total: 10, frames_rejected: 0 })))
      .toBe("Beamline average of 10 images");
    expect(exposureKindLabel(makeExposure({ frames_total: 1, frames_rejected: 0 })))
      .toBe("Beamline average of 1 image");
  });

  it("reads as a single image while the counts are unknown", () => {
    expect(exposureKindLabel(makeExposure({ frames_total: null, frames_rejected: null }))).toBe("Single");
  });
});

describe("exposureTimeLabel", () => {
  it("shows the exposure time of one image", () => {
    expect(exposureTimeLabel(makeExposure({ exposure_time: 1 }))).toBe("1.0 s");
    expect(exposureTimeLabel(makeExposure())).toBe("—");
  });

  it("multiplies the time of one image by the images the beamline kept", () => {
    expect(exposureTimeLabel(makeExposure({ exposure_time: 1, frames_total: 10, frames_rejected: 1 })))
      .toBe("9 × 1.0 s");
    expect(exposureTimeLabel(makeExposure({ exposure_time: 0.5, frames_total: 10, frames_rejected: 0 })))
      .toBe("10 × 0.5 s");
  });

  it("shows a dash when an average has no exposure time", () => {
    expect(exposureTimeLabel(makeExposure({ frames_total: 10, frames_rejected: 0 }))).toBe("—");
  });
});
