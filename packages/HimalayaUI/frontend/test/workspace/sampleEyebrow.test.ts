import { describe, expect, it } from "vitest";
import { sampleEyebrow } from "../../src/workspace/sampleEyebrow";

describe("sampleEyebrow", () => {
  it("names the setup and the sample's place in the list", () => {
    expect(sampleEyebrow("1.7 m", 2, 13)).toBe("1.7 m · sample 3 of 13");
    expect(sampleEyebrow(undefined, 0, 4)).toBe("Sample 1 of 4");
    expect(sampleEyebrow("1.7 m", -1, 4)).toBe("1.7 m");
  });
});
