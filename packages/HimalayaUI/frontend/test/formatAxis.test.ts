/**
 * Tests for the shared plain-decimal axis formatter (issue #80). Pins the
 * SAXS-friendly tick vocabulary: scientific notation only at the extremes,
 * plain decimals in the readable middle band, and absolute-value branching
 * so negatives get the same tier as their positive twin. Both `TraceViewer`
 * (Index page) and `MultiTracePlot` (Compare page) consume this — keeping
 * the cases here means a single source of truth for what a tick looks like.
 */
import { describe, it, expect } from "vitest";
import { formatAxis } from "../src/lib/plot/formatAxis";

describe("formatAxis", () => {
	it("returns '0' for exact zero", () => {
		expect(formatAxis(0)).toBe("0");
	});

	it("uses scientific notation for very small magnitudes (< 1e-3)", () => {
		// 0.0001 = 1e-4 → scientific notation branch (toExponential(0) → "1e-4").
		expect(formatAxis(0.0001)).toBe("1e-4");
	});

	it("uses 3 dp for the 0.001 ≤ |d| < 0.01 band", () => {
		expect(formatAxis(0.005)).toBe("0.005");
	});

	it("uses 2 dp for the 0.01 ≤ |d| < 1 band (small)", () => {
		expect(formatAxis(0.05)).toBe("0.05");
	});

	it("uses 2 dp for the 0.01 ≤ |d| < 1 band (mid)", () => {
		expect(formatAxis(0.5)).toBe("0.50");
	});

	it("uses 1 dp for the 1 ≤ |d| < 10 band", () => {
		expect(formatAxis(5)).toBe("5.0");
	});

	it("uses 0 dp for the 10 ≤ |d| < 100 band", () => {
		expect(formatAxis(50)).toBe("50");
	});

	it("uses 0 dp for the 100 ≤ |d| < 1e4 band", () => {
		expect(formatAxis(500)).toBe("500");
	});

	it("uses scientific notation for very large magnitudes (≥ 1e4)", () => {
		// 1e5 → "1e5" after stripping "+0" / "+".
		expect(formatAxis(1e5)).toBe("1e5");
	});

	it("preserves the sign on negative inputs (absolute-value path)", () => {
		expect(formatAxis(-0.05)).toBe("-0.05");
	});

	it("preserves the sign on negative scientific-notation inputs", () => {
		expect(formatAxis(-0.0001)).toBe("-1e-4");
	});
});
