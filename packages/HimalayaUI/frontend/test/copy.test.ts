import { describe, it, expect } from "vitest";
import { sanitizeDashes } from "../src/lib/copy";

describe("sanitizeDashes", () => {
  it("folds a spaced em dash from upstream data into the middot separator", () => {
    expect(sanitizeDashes("SSRL April 2026 — 1p7m")).toBe("SSRL April 2026 · 1p7m");
  });

  it("folds en dashes too", () => {
    expect(sanitizeDashes("10–20 keV")).toBe("10 · 20 keV");
  });

  it("leaves clean copy (hyphens, middots) untouched", () => {
    expect(sanitizeDashes("C1 · SSRL · exp 5")).toBe("C1 · SSRL · exp 5");
    expect(sanitizeDashes("low-amplitude")).toBe("low-amplitude");
  });

  it("is idempotent", () => {
    const once = sanitizeDashes("a — b — c");
    expect(sanitizeDashes(once)).toBe(once);
  });

  it("emits no em dash on any folded input", () => {
    expect(sanitizeDashes("x—y–z")).not.toMatch(/[—–]/);
  });
});
