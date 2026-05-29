import { describe, it, expect } from "vitest";
import { scanContent, hashViolation, diffBaseline } from "../scripts/check-design.mjs";

// Helper: collect the set of rule ids a single line trips (outside ui/).
function rulesFor(line: string): string[] {
  return scanContent("components/Probe.tsx", line).map((v) => v.rule);
}

describe("check-design guard — ban rules (spec §4)", () => {
  // Rule #1 — arbitrary text size/color
  it("flags text-[10px]", () => {
    expect(rulesFor('<span className="text-[10px] mt-1" />')).toContain("no-arbitrary-text");
  });
  it("passes the on-scale text-base", () => {
    expect(rulesFor('<span className="text-base font-medium" />')).not.toContain("no-arbitrary-text");
  });

  // Rule #2 — arbitrary radius
  it("flags rounded-[3px]", () => {
    expect(rulesFor('<div className="rounded-[3px]" />')).toContain("no-arbitrary-radius");
  });
  it("passes rounded-md", () => {
    expect(rulesFor('<div className="rounded-md bg-plate" />')).not.toContain("no-arbitrary-radius");
  });

  // Rule #4 — side stripe
  it("flags border-l-4", () => {
    expect(rulesFor('<div className="border-l-4 border-error" />')).toContain("no-side-stripe");
  });

  // Rule #3 — raw color in a color utility, shadow-[…] excluded
  it("flags bg-[oklch(...)] in a non-allowlisted file", () => {
    expect(rulesFor('<div className="bg-[oklch(0.05_0_0/0.65)]" />')).toContain("no-raw-color-utility");
  });
  it("passes a shadow-[…rgba…] Plate-Lift literal (rule #3 excludes shadow-)", () => {
    expect(rulesFor('<div className="shadow-[0_8px_26px_-10px_rgba(60,52,40,.34)]" />')).not.toContain("no-raw-color-utility");
  });
  it("allowlists bg-[oklch] inside an allowlisted file (rule #3 only)", () => {
    const v = scanContent("components/DetectorImage.tsx", '<div className="bg-[oklch(0.15_0.01_55)]" />');
    expect(v.map((x) => x.rule)).not.toContain("no-raw-color-utility");
  });

  // Rule #5 — raw color literal anywhere; var(--color-*) stripped first
  it("passes a var-only color-mix string literal", () => {
    const line = '  background: "color-mix(in oklab, var(--color-accent) 9%, transparent)",';
    expect(rulesFor(line)).not.toContain("no-raw-color-literal");
  });
  it("passes a computed style={{ color }} identifier", () => {
    expect(rulesFor("        style={{ color }}")).not.toContain("no-raw-color-literal");
  });
  it("passes a template-literal interpolated color", () => {
    expect(rulesFor("        style={{ background: `color-mix(in oklab, ${x} 20%, transparent)` }}")).not.toContain("no-raw-color-literal");
  });
  it("flags a raw oklch scrim string literal", () => {
    expect(rulesFor('        background: "oklch(0.05 0 0 / 0.65)",')).toContain("no-raw-color-literal");
  });
  it("flags a bare var() with NO raw literal as clean (FocusReflectionsTable/SeriesFolioCard cases)", () => {
    expect(rulesFor('        style={{ background: "var(--color-success)" }}')).not.toContain("no-raw-color-literal");
  });
  it("passes a PR-number cross-reference in a comment (not a color)", () => {
    expect(rulesFor("  // see #183 and #221 for the rationale")).not.toContain("no-raw-color-literal");
  });
  it("flags a quoted raw hex color literal", () => {
    expect(rulesFor('        style={{ color: "#3c3428" }}')).toContain("no-raw-color-literal");
  });

  // Scope exclusion — src/components/ui/** is never scanned
  it("excludes src/components/ui/** entirely", () => {
    expect(scanContent("components/ui/Toast.tsx", '<div className="border-l-4 text-[10px]" />')).toHaveLength(0);
  });
});

describe("check-design guard — baseline diff", () => {
  it("returns clean (no new) when every violation hash is in the baseline", () => {
    const violations = scanContent("components/Probe.tsx", '<div className="text-[10px]" />');
    const baseline = { hashes: violations.map((v) => hashViolation(v.rule, v.text)) };
    const { notInBaseline, grew } = diffBaseline(violations, baseline);
    expect(notInBaseline).toHaveLength(0);
    expect(grew).toBe(false);
  });

  it("flags a not-in-baseline violation (CI would exit 2)", () => {
    const violations = scanContent("components/Probe.tsx", '<div className="text-[10px]" />');
    const { notInBaseline, grew } = diffBaseline(violations, { hashes: [] });
    expect(notInBaseline).toHaveLength(1);
    expect(grew).toBe(true);
  });

  it("hash is stable across whitespace/indent reflow (move-between-files is a no-op)", () => {
    expect(hashViolation("no-arbitrary-text", "text-[10px]")).toBe(
      hashViolation("no-arbitrary-text", "   text-[10px]   "),
    );
  });
});
