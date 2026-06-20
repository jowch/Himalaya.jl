import { describe, it, expect } from "vitest";
import { matchSample } from "../src/lib/matchSample";
import type { LoadSample } from "../src/api";

const s: LoadSample = {
  sample_id: 1, name: "HA85 (S01P15)", slot_index: 15,
  grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null,
  exposures: [{ id: 9, filename: "HA_0612_001.tif", horizontal_position: 1, timestamp: null }],
};

describe("matchSample", () => {
  it("empty query matches everything", () => { expect(matchSample(s, "")).toBe(true); });
  it("substring matches name (case-insensitive)", () => { expect(matchSample(s, "ha85")).toBe(true); });
  it("substring matches an exposure filename", () => { expect(matchSample(s, "0612")).toBe(true); });
  it("glob * matches a prefix", () => { expect(matchSample(s, "ha8*")).toBe(true); });
  it("glob ? matches one char", () => { expect(matchSample({ ...s, name: "JC C04" }, "JC C0?")).toBe(true); });
  it("non-matching glob fails", () => { expect(matchSample(s, "ZZ*")).toBe(false); });
});
