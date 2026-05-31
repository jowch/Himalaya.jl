import { describe, it, expect } from "vitest";
import { scanContent, scanLegacyImports } from "../scripts/check-design.mjs";

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
  it("flags a bare var() with NO raw literal as clean (SeriesFolioCard case)", () => {
    expect(rulesFor('        style={{ background: "var(--color-success)" }}')).not.toContain("no-raw-color-literal");
  });
  it("passes a PR-number cross-reference in a comment (not a color)", () => {
    expect(rulesFor("  // see #183 and #221 for the rationale")).not.toContain("no-raw-color-literal");
  });
  it("flags a quoted raw hex color literal", () => {
    expect(rulesFor('        style={{ color: "#3c3428" }}')).toContain("no-raw-color-literal");
  });

  // Rule #5 — shares rule #3's color-authoring allowlist
  it("does NOT flag a raw oklch literal in an allowlisted color-authoring file (rule #5)", () => {
    const line = '  Pn3m: { hue: "oklch(0.62 0.18 28)" },';
    expect(scanContent("phases.ts", line).map((v) => v.rule)).not.toContain("no-raw-color-literal");
  });
  it("DOES flag the same raw oklch literal in a normal consumer (rule #5)", () => {
    const line = '  Pn3m: { hue: "oklch(0.62 0.18 28)" },';
    expect(scanContent("components/Probe.tsx", line).map((v) => v.rule)).toContain("no-raw-color-literal");
  });

  // Rule #5 — shadow-[…rgba…] elevation tokens are stripped first (consistent with rule #3,
  // which excludes shadow-): the rgba lives inside an elevation token, not a color role.
  it("passes a shadow-[…rgba…] Plate-Lift literal under rule #5 too", () => {
    const line = '      className="rounded-md shadow-[0_8px_26px_-10px_rgba(60,52,40,.34)]"';
    expect(rulesFor(line)).not.toContain("no-raw-color-literal");
  });

  // Scope exclusion — src/components/ui/** is never scanned
  it("excludes src/components/ui/** entirely", () => {
    expect(scanContent("components/ui/Toast.tsx", '<div className="border-l-4 text-[10px]" />')).toHaveLength(0);
  });
});

describe("scanLegacyImports", () => {
  it("flags a print/ file importing old components (one level up to src)", () => {
    const v = scanLegacyImports("print/App.tsx", `import { X } from "../components/X";`);
    expect(v).toEqual([
      { rule: "no-legacy-import", file: "print/App.tsx", line: 1, text: "../components/X" },
    ]);
  });

  it("flags a print/ui file importing old components (two levels up to src)", () => {
    const v = scanLegacyImports("print/ui/Button.tsx", `import { X } from "../../components/X";`);
    expect(v.map((h) => h.text)).toEqual(["../../components/X"]);
  });

  it("flags importing old pages", () => {
    const v = scanLegacyImports("print/pages/P.tsx", `import P from "../../pages/Old";`);
    expect(v.map((h) => h.rule)).toEqual(["no-legacy-import"]);
  });

  it("flags importing the OLD ui primitives (under components/)", () => {
    const v = scanLegacyImports("print/ui/Card.tsx", `import { Dot } from "../../components/ui/Dot";`);
    expect(v).toHaveLength(1);
  });

  it("allows importing the shared core", () => {
    const src = [
      `import { api } from "../api";`,
      `import { phaseColor } from "../phases";`,
      `import { foo } from "../lib/queue/types";`,
      `import { useFocusTrap } from "../hooks/useFocusTrap";`,
    ].join("\n");
    expect(scanLegacyImports("print/App.tsx", src)).toEqual([]);
  });

  it("allows print-internal imports (its own components/pages)", () => {
    const v = scanLegacyImports("print/App.tsx", `import { Page } from "./pages/ContactSheet";\nimport { C } from "./components/C";`);
    expect(v).toEqual([]);
  });

  it("ignores non-print files entirely", () => {
    expect(scanLegacyImports("components/Foo.tsx", `import x from "../pages/Bar";`)).toEqual([]);
  });
});

describe("isExcluded via scanContent — print/ui authoring", () => {
  it("excludes src/print/ui from appearance rules", () => {
    expect(scanContent("print/ui/PhaseChip.tsx", `<span className="text-[11px]" />`)).toEqual([]);
  });
  it("still enforces appearance rules on print/ consumers", () => {
    const v = scanContent("print/components/Foo.tsx", `<span className="text-[11px]" />`);
    expect(v.map((h) => h.rule)).toContain("no-arbitrary-text");
  });
});
