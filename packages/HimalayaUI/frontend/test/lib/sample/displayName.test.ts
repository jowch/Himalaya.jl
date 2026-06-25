import { describe, it, expect } from "vitest";
import { sampleDisplayName } from "../../../src/lib/sample/displayName";

describe("sampleDisplayName", () => {
  it("returns name when set", () => {
    expect(sampleDisplayName({ id: 1, name: "HA85 (S01P15)" })).toBe("HA85 (S01P15)");
  });
  it("falls back to `Sample id` when name is empty string (uses ||, not ??)", () => {
    expect(sampleDisplayName({ id: 7, name: "" })).toBe("Sample 7");
  });
});
