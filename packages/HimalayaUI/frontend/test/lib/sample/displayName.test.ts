import { describe, it, expect } from "vitest";
import { sampleDisplayName } from "../../../src/lib/sample/displayName";

describe("sampleDisplayName", () => {
  it("prefers display_name when set", () => {
    expect(sampleDisplayName({ id: 1, name: "JC001", display_name: "DOPC" })).toBe("DOPC");
  });
  it("falls back to name when display_name is null", () => {
    expect(sampleDisplayName({ id: 1, name: "JC001", display_name: null })).toBe("JC001");
  });
  it("falls back to name when display_name is empty string (uses ||, not ??)", () => {
    expect(sampleDisplayName({ id: 1, name: "JC001", display_name: "" })).toBe("JC001");
  });
  it("falls back to `Sample #id` when both are null/empty", () => {
    expect(sampleDisplayName({ id: 7, name: null, display_name: null })).toBe("Sample #7");
    expect(sampleDisplayName({ id: 7, name: "", display_name: "" })).toBe("Sample #7");
  });
});
