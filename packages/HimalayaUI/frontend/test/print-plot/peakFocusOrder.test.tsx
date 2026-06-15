import { describe, it, expect } from "vitest";
import { nextFocusPeakId } from "../../src/print/plot/peakFocusOrder";

describe("nextFocusPeakId", () => {
  const peaks = [
    { id: 10, q: 0.05 },
    { id: 20, q: 0.1 },
    { id: 30, q: 0.2 },
  ];

  it("prefers the next peak by q-order", () => {
    expect(nextFocusPeakId(peaks, 10)).toBe(20);
    expect(nextFocusPeakId(peaks, 20)).toBe(30);
  });

  it("falls back to the previous peak when the removed one is the last (highest q)", () => {
    expect(nextFocusPeakId(peaks, 30)).toBe(20);
  });

  it("returns null when no peak would survive", () => {
    expect(nextFocusPeakId([{ id: 5, q: 0.1 }], 5)).toBeNull();
  });

  it("returns null when the removed id is not present", () => {
    expect(nextFocusPeakId(peaks, 999)).toBeNull();
  });

  it("orders by q, not by array insertion order (manual peaks land out of order)", () => {
    // id 30 (q=0.2) sits AFTER id 20 (q=0.1) by q despite being earlier in the
    // array → removing 30 should hand focus back to 20, the q-previous peak.
    const outOfOrder = [
      { id: 30, q: 0.2 },
      { id: 10, q: 0.05 },
      { id: 20, q: 0.1 },
    ];
    expect(nextFocusPeakId(outOfOrder, 10)).toBe(20); // q-next of 0.05 is 0.1
    expect(nextFocusPeakId(outOfOrder, 30)).toBe(20); // 0.2 is highest → q-prev 0.1
  });
});
