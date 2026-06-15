import { describe, it, expect } from "vitest";
import { reorder, dropEdgeFor } from "../../src/print/components/useDragReorder";
describe("dropEdgeFor", () => {
  it("is null when not dragging", () => expect(dropEdgeFor(null, 2, 2)).toBeNull());
  it("is null off the over-row", () => expect(dropEdgeFor(0, 2, 3)).toBeNull());
  it("is null on the dragged row itself", () => expect(dropEdgeFor(2, 2, 2)).toBeNull());
  it("shows TOP when dragging up (from below the target)", () => expect(dropEdgeFor(3, 1, 1)).toBe("top"));
  it("shows BOTTOM when dragging down (from above the target)", () => expect(dropEdgeFor(1, 3, 3)).toBe("bottom"));
});
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
