import { describe, expect, it } from "vitest";
import { describeRibbon, reviewSummary, ribbonAt, ribbonLayout, toReview } from "../../src/lib/collectionRibbon";
import type { CollectionTimeline } from "../../src/api";

const TL: CollectionTimeline = {
  start: "2026-04-25T10:00:00",
  samples: [[0, 300, 1], [420, 900, 0], [4500, 5400, 0]],
  beam_off: [5400],
};

describe("ribbonLayout", () => {
  it("joins samples measured within ten minutes of each other into one stretch", () => {
    const { spans, end } = ribbonLayout(TL);
    expect(end).toBe(5400);
    expect(spans.map((s) => [s.from, s.to, s.reviewed])).toEqual([[0, 420, true], [420, 900, false], [4500, 5400, false]]);
  });

  it("leaves a gap of more than ten minutes empty", () => {
    const { gaps } = ribbonLayout(TL);
    expect(gaps).toEqual([{ from: 900, to: 4500 }]);
  });
});

describe("ribbonAt", () => {
  it("names the sample under the pointer and whether it is reviewed", () => {
    expect(ribbonAt(TL, 100)).toEqual({ kind: "sample", from: 0, to: 300, reviewed: true, beamOff: 0 });
  });

  it("counts images taken with the beam off in that sample", () => {
    expect(ribbonAt(TL, 5000)).toMatchObject({ kind: "sample", beamOff: 1 });
  });

  it("reports how long a gap lasted", () => {
    expect(ribbonAt(TL, 2000)).toEqual({ kind: "gap", from: 900, to: 4500 });
  });
});

describe("describeRibbon", () => {
  it("gives the start and the length", () => {
    expect(describeRibbon(TL)).toEqual({ started: "Sat 25 Apr 10:00", length: "1 h 30 min" });
  });
});

describe("toReview", () => {
  it("gives how many samples are left to review", () => {
    expect(toReview({ reviewed: 1, total: 3 })).toBe("2 to review");
  });

  it("says when every sample is reviewed", () => {
    expect(toReview({ reviewed: 3, total: 3 })).toBe("All reviewed");
    expect(reviewSummary({ reviewed: 1, total: 1 })).toBe("1 of 1 sample reviewed");
  });
});
