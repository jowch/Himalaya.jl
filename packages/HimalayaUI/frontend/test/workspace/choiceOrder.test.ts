import { describe, expect, it } from "vitest";
import type { SampleLabel } from "../../src/api";
import { backgroundOrder, compareOrder } from "../../src/workspace/choiceOrder";

const sample = (id: number, name: string): { id: number; name: string } => ({ id, name });
const buffer = (sampleId: number, state: SampleLabel["state"] = "confirmed"): SampleLabel => ({ id: sampleId, sample_id: sampleId, label: "Buffer", state, source: "user" });
const names = (xs: readonly { name: string }[]): string[] => xs.map((x) => x.name);

describe("compareOrder", () => {
  const samples = [sample(1, "X_10A"), sample(2, "X_02A"), sample(3, "buf_2"), sample(4, "X_01B"), sample(5, "buf_10")];

  it("lists samples in natural name order when nothing was used and nothing is a buffer", () => {
    expect(names(compareOrder(samples.filter((s) => s.id !== 3 && s.id !== 5), [], []))).toEqual(["X_01B", "X_02A", "X_10A"]);
  });

  it("puts recently compared samples first, most recent first, then buffers, then the rest by name", () => {
    expect(names(compareOrder(samples, [buffer(3), buffer(5, "suggested")], ["1", "4"]))).toEqual(["X_10A", "X_01B", "buf_2", "buf_10", "X_02A"]);
  });

  it("orders buffers of the same standing by name", () => {
    expect(names(compareOrder(samples, [buffer(5), buffer(3)], []))).toEqual(["buf_2", "buf_10", "X_01B", "X_02A", "X_10A"]);
  });
});

describe("backgroundOrder", () => {
  const samples = [sample(1, "S_10"), sample(2, "S_02"), sample(3, "B_run"), sample(4, "B_setup"), sample(5, "S_01"), sample(6, "open")];

  it("lists the run's buffers, then the setup's, then every other sample by name, leaving out the open sample and its background", () => {
    const out = backgroundOrder({ samples, open: 6, current: 5, runBuffers: [3], setupBuffers: [4, 3], recent: [] });
    expect(names(out)).toEqual(["B_run", "B_setup", "S_02", "S_10"]);
  });

  it("puts backgrounds chosen recently ahead of the run's buffers, most recent first", () => {
    const out = backgroundOrder({ samples, open: 6, current: null, runBuffers: [3], setupBuffers: [4], recent: ["5", "1", "6"] });
    expect(names(out)).toEqual(["S_01", "S_10", "B_run", "B_setup", "S_02"]);
  });
});
