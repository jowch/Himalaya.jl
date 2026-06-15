// test/newSeriesNav.test.ts
import { describe, it, expect, vi } from "vitest";
import { navigateToNewSeries, readNewSeriesSeed } from "../src/lib/series/newSeriesNav";

describe("navigateToNewSeries", () => {
  it("calls navigate with /series/new and the seed as location state", () => {
    const navigate = vi.fn();
    navigateToNewSeries(new Set([10, 11, 12]), navigate);
    expect(navigate).toHaveBeenCalledWith("/series/new", {
      state: { seedSampleIds: expect.arrayContaining([10, 11, 12]) },
    });
    // Deterministic: always an array (not a Set) so it round-trips through history.
    const call = navigate.mock.calls[0]!;
    expect(Array.isArray((call[1] as { state: { seedSampleIds: unknown } }).state.seedSampleIds)).toBe(true);
  });
});

describe("readNewSeriesSeed", () => {
  it("returns the seed array when location state is well-formed", () => {
    const loc = { state: { seedSampleIds: [10, 11] } } as unknown as Location;
    expect(readNewSeriesSeed(loc)).toEqual([10, 11]);
  });

  it("returns null when location has no state", () => {
    expect(readNewSeriesSeed({ state: null } as unknown as Location)).toBeNull();
  });

  it("returns null when seedSampleIds is not an array", () => {
    expect(readNewSeriesSeed({ state: { seedSampleIds: "bad" } } as unknown as Location)).toBeNull();
  });

  it("returns null when seedSampleIds is missing", () => {
    expect(readNewSeriesSeed({ state: {} } as unknown as Location)).toBeNull();
  });
});
