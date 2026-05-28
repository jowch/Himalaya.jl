import { describe, it, expect } from "vitest";
import { parseSortKey } from "../src/lib/scoping/parseSortKey";

describe("parseSortKey", () => {
  it("parses a plain number", () => {
    expect(parseSortKey("2.5")).toBe(2.5);
    expect(parseSortKey("0")).toBe(0);
  });
  it("parses the divisor of an a : b ratio (low→high by the second term)", () => {
    expect(parseSortKey("1 : 0")).toBe(0);
    expect(parseSortKey("1 : 0.25")).toBe(0.25);
    expect(parseSortKey("1:4")).toBe(4);
  });
  it("parses a trailing number embedded in text", () => {
    expect(parseSortKey("dose 30mM")).toBe(30);
  });
  it("returns null for an unparseable value", () => {
    expect(parseSortKey("DOPC")).toBeNull();
    expect(parseSortKey("")).toBeNull();
  });
});
