import { describe, expect, it } from "vitest";
import { rankBridgeMates } from "../../src/lib/bridges";

const s = (id: number, name: string, setup_id: number | null): { id: number; name: string; setup_id: number | null } => ({ id, name, setup_id });

describe("rankBridgeMates", () => {
  it("puts samples from other setups first, the same name before a similar one, and keeps the list order otherwise", () => {
    const here = s(1, "MO_H2O", 10);
    const others = [s(2, "MO_H2O_37C", 10), s(3, "DOPC", 20), s(4, "MO_PEG5", 20), s(5, "mo-h2o", 20), s(6, "MO_H2O_2", 20), s(7, "MO_H2O", 10)];
    expect(rankBridgeMates(here, others).map((x) => x.id)).toEqual([5, 6, 4, 3, 7, 2]);
  });

  it("ranks by name alone when the sample has no setup", () => {
    const others = [s(2, "B", null), s(3, "A", null)];
    expect(rankBridgeMates(s(1, "A", null), others).map((x) => x.id)).toEqual([3, 2]);
  });
});
