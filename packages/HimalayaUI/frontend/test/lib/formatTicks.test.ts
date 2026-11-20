import { describe, it, expect } from "vitest";
import { formatAxis, formatTicks } from "../../src/lib/plot/formatAxis";

describe("formatTicks", () => {
  it("keeps the plain labels when they are already distinct", () => {
    const values = [0.01, 0.1, 1, 10];
    expect(formatTicks(values)).toEqual(values.map(formatAxis));
  });

  it("drops the decimal when every tick is a whole number", () => {
    expect(formatTicks([5, 10, 15])).toEqual(["5", "10", "15"]);
  });

  it("adds digits until close ticks get distinct labels", () => {
    expect(formatTicks([0.1, 0.105, 0.11, 0.115])).toEqual(["0.100", "0.105", "0.110", "0.115"]);
  });
});
