import { describe, it, expect } from "vitest";
import type { SampleBridge } from "../../src/api";
import { applyBridge, applyUnbridge, bridgedWith, bridgeCounts } from "../../src/lib/bridges";

const b = (sample_id: number, group_id: number): SampleBridge => ({ sample_id, group_id });

function groups(list: readonly SampleBridge[]): number[][] {
  const by = new Map<number, number[]>();
  for (const r of list) by.set(r.group_id, [...(by.get(r.group_id) ?? []), r.sample_id]);
  return [...by.values()].map((g) => [...g].sort((x, y) => x - y)).sort((x, y) => x[0]! - y[0]!);
}

describe("bridges", () => {
  it("lists the other samples in a sample's group", () => {
    const list = [b(1, 7), b(2, 7), b(3, 7), b(4, 8), b(5, 8)];
    expect(bridgedWith(list, 2)).toEqual([1, 3]);
    expect(bridgedWith(list, 9)).toEqual([]);
  });

  it("counts the other members for each bridged sample", () => {
    const counts = bridgeCounts([b(1, 7), b(2, 7), b(3, 7), b(4, 8), b(5, 8)]);
    expect(counts.get(1)).toBe(2);
    expect(counts.get(4)).toBe(1);
    expect(counts.has(6)).toBe(false);
  });

  it("bridging two loose samples makes a new group", () => {
    expect(groups(applyBridge([b(4, 8), b(5, 8)], 1, 2))).toEqual([[1, 2], [4, 5]]);
  });

  it("bridging joins an existing group from either side", () => {
    expect(groups(applyBridge([b(1, 7), b(2, 7)], 3, 1))).toEqual([[1, 2, 3]]);
    expect(groups(applyBridge([b(1, 7), b(2, 7)], 1, 3))).toEqual([[1, 2, 3]]);
  });

  it("bridging across two groups combines them", () => {
    expect(groups(applyBridge([b(1, 7), b(2, 7), b(3, 8), b(4, 8)], 2, 4))).toEqual([[1, 2, 3, 4]]);
  });

  it("bridging a sample to itself or within its group changes nothing", () => {
    const list = [b(1, 7), b(2, 7)];
    expect(applyBridge(list, 1, 1)).toEqual(list);
    expect(applyBridge(list, 1, 2)).toEqual(list);
  });

  it("unbridging takes the sample out, and a group of one dissolves", () => {
    expect(groups(applyUnbridge([b(1, 7), b(2, 7), b(3, 7)], 3))).toEqual([[1, 2]]);
    expect(applyUnbridge([b(1, 7), b(2, 7), b(4, 8), b(5, 8)], 1)).toEqual([b(4, 8), b(5, 8)]);
    expect(applyUnbridge([b(1, 7), b(2, 7)], 9)).toEqual([b(1, 7), b(2, 7)]);
  });
});
