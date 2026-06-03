import { describe, it, expect } from "vitest";
import { reorder } from "../../src/print/components/useDragReorder";
describe("reorder", () => {
  it("moves an item down", () => expect(reorder(["a", "b", "c", "d"], 0, 2)).toEqual(["b", "c", "a", "d"]));
  it("moves an item up", () => expect(reorder(["a", "b", "c", "d"], 3, 1)).toEqual(["a", "d", "b", "c"]));
  it("is a no-op for same index", () => expect(reorder(["a", "b", "c"], 1, 1)).toEqual(["a", "b", "c"]));
  it("returns a copy (no mutation)", () => {
    const src = ["a", "b"];
    const out = reorder(src, 0, 1);
    expect(out).not.toBe(src);
    expect(src).toEqual(["a", "b"]);
  });
  it("ignores out-of-range", () => expect(reorder(["a", "b"], 0, 5)).toEqual(["a", "b"]));
});
