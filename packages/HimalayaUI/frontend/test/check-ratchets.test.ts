import { describe, it, expect } from "vitest";
import {
  countAppearance,
  countTsComments,
  countJuliaComments,
  compareToBaseline,
} from "../scripts/check-ratchets.mjs";

describe("countAppearance", () => {
  it("counts each className that carries appearance tokens", () => {
    const src = [
      '<span className="text-caption text-ink-soft font-mono" />',
      '<div className="mt-2 flex gap-1" />',
      '<p className={`text-meta ${x}`} />',
    ].join("\n");
    expect(countAppearance("print/pages/X.tsx", src)).toBe(2);
  });

  it("ignores primitives, render layers and stories", () => {
    const src = '<span className="text-caption text-ink" />';
    expect(countAppearance("print/ui/Text.tsx", src)).toBe(0);
    expect(countAppearance("print/plot/TracePlot.tsx", src)).toBe(0);
    expect(countAppearance("print/pages/X.stories.tsx", src)).toBe(0);
  });

  it("treats text alignment and wrapping as placement", () => {
    expect(countAppearance("print/pages/X.tsx", '<p className="text-center text-nowrap" />')).toBe(0);
  });
});

describe("countTsComments", () => {
  it("counts line, block and JSX comment lines", () => {
    const src = [
      "// one",
      "const a = 1;",
      "/**",
      " * two",
      " */",
      "{/* three */}",
    ].join("\n");
    expect(countTsComments(src)).toBe(5);
  });

  it("exempts tool directives and physics rationale", () => {
    const src = [
      "// eslint-disable-next-line no-console",
      "// @ts-expect-error upstream types",
      "// physics: q = 4π sin(θ)/λ",
    ].join("\n");
    expect(countTsComments(src)).toBe(0);
  });
});

describe("countJuliaComments", () => {
  it("counts # lines and top-level docstrings", () => {
    const src = [
      '"""',
      "    f(x)",
      "Doc.",
      '"""',
      "function f(x)",
      "    # inline note",
      "    x",
      "end",
    ].join("\n");
    expect(countJuliaComments(src)).toBe(5);
  });

  it("does not count indented SQL strings or physics lines", () => {
    const src = [
      "rows = DBInterface.execute(db,",
      '    """SELECT id',
      "       FROM t",
      '    """)',
      "# physics: Q0f31 detector counts are log-normalised",
    ].join("\n");
    expect(countJuliaComments(src)).toBe(0);
  });
});

describe("compareToBaseline", () => {
  const baseline = { appearance: { "a.tsx": 3 }, comments: { "a.tsx": 10 } };

  it("passes when counts match", () => {
    const r = compareToBaseline(baseline, { appearance: { "a.tsx": 3 }, comments: { "a.tsx": 10 } });
    expect(r.increased).toEqual([]);
    expect(r.decreased).toEqual([]);
  });

  it("reports a file whose count grew, and a new file with any count", () => {
    const r = compareToBaseline(baseline, {
      appearance: { "a.tsx": 4, "b.tsx": 1 },
      comments: { "a.tsx": 10 },
    });
    expect(r.increased.map((x) => `${x.kind}:${x.file}`)).toEqual(["appearance:a.tsx", "appearance:b.tsx"]);
  });

  it("reports a file whose count fell so the baseline can be tightened", () => {
    const r = compareToBaseline(baseline, { appearance: { "a.tsx": 1 }, comments: {} });
    expect(r.decreased.map((x) => `${x.kind}:${x.file}`)).toEqual(["appearance:a.tsx", "comments:a.tsx"]);
  });
});
