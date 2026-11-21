import { describe, it, expect } from "vitest";
import { formatAxis, formatTicks } from "../../src/lib/plot/formatAxis";

describe("formatTicks", () => {
  it("keeps the plain labels when they are distinct and exact", () => {
    const values = [0.01, 0.1, 1, 10];
    expect(formatTicks(values)).toEqual(values.map(formatAxis));
  });

  it("drops the decimal when every tick is a whole number", () => {
    expect(formatTicks([5, 10, 15])).toEqual(["5", "10", "15"]);
    expect(formatTicks([5])).toEqual(["5"]);
  });

  it("keeps decades from 1e4 up in exponent form", () => {
    expect(formatTicks([1000, 10000, 100000])).toEqual(["1000", "1e4", "1e5"]);
  });

  it("adds only the digits close ticks need", () => {
    expect(formatTicks([0.1, 0.105, 0.11, 0.115])).toEqual(["0.100", "0.105", "0.110", "0.115"]);
    expect(formatTicks([12.5, 15, 17.5])).toEqual(["12.5", "15.0", "17.5"]);
  });

  it("never rounds a large or small tick to a neighbouring value", () => {
    expect(formatTicks([5000, 10000, 15000])).toEqual(["5000", "1e4", "1.5e4"]);
    expect(formatTicks([20000, 20200, 20400])).toEqual(["2e4", "2.02e4", "2.04e4"]);
    expect(formatTicks([1.2e-4, 1.4e-4, 1.6e-4])).toEqual(["1.2e-4", "1.4e-4", "1.6e-4"]);
  });
});
