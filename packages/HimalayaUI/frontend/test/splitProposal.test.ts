import { describe, it, expect } from "vitest";
import { splitProposal, humanizeKey } from "../src/lib/scoping/splitProposal";
import type { OrderingProposal } from "../src/lib/scoping/proposeOrdering";

const proposal = (rows: OrderingProposal["rows"], key: string | undefined): OrderingProposal =>
  ({ orderingKey: key, rows });

describe("splitProposal", () => {
  it("members carry a parsed value; loose matches lack the key (value empty)", () => {
    const p = proposal([
      { sampleId: 10, sampleName: "A", value: "1:1", flagged: false, include: true },
      { sampleId: 11, sampleName: "B", value: "", flagged: true, include: true },
    ], "ratio");
    const { members, looseMatches } = splitProposal(p);
    expect(members.map((r) => r.sampleId)).toEqual([10]);
    expect(looseMatches.map((r) => r.sampleId)).toEqual([11]);
  });
  it("loose matches start excluded (include=false)", () => {
    const p = proposal([
      { sampleId: 11, sampleName: "B", value: "", flagged: true, include: true },
    ], "ratio");
    expect(splitProposal(p).looseMatches[0].include).toBe(false);
  });
  it("a cold corpus (no orderingKey) puts every row in looseMatches", () => {
    const p = proposal([
      { sampleId: 10, sampleName: "A", value: "", flagged: true, include: true },
    ], undefined);
    const { members, looseMatches } = splitProposal(p);
    expect(members).toEqual([]);
    expect(looseMatches).toHaveLength(1);
  });
});

describe("humanizeKey", () => {
  it("turns snake/kebab keys into spaced words", () => {
    expect(humanizeKey("ll37_lipid_ratio")).toBe("ll37 lipid ratio");
    expect(humanizeKey("ratio")).toBe("ratio");
    expect(humanizeKey(undefined)).toBe("—");
  });
});
