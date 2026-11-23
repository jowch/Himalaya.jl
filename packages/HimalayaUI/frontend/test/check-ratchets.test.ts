import { describe, it, expect } from "vitest";
import {
  countAppearance,
  countTsComments,
  countJuliaComments,
  compareToBaseline,
  canUpdate,
  TS_ROOTS,
  JULIA_ROOTS,
} from "../scripts/check-ratchets.mjs";

describe("countAppearance", () => {
  it("counts every appearance token, not just styled classNames", () => {
    const src = [
      '<span className="text-caption text-ink-soft font-mono" />',
      '<div className="mt-2 flex gap-1" />',
      '<p className={`text-meta ${x}`} />',
    ].join("\n");
    expect(countAppearance("pages/X.tsx", src)).toBe(4);
  });

  it("reads every literal inside cx() calls, ternaries and template expressions", () => {
    expect(countAppearance("pages/X.tsx", '<i className={cx("mt-2", on && "text-ink-soft font-mono")} />')).toBe(2);
    expect(countAppearance("pages/X.tsx", '<i className={x ? "mt-1" : "bg-paper"} />')).toBe(1);
    expect(countAppearance("pages/X.tsx", '<i className={`mt-1 ${on ? "uppercase" : ""}`} />')).toBe(1);
  });

  it("grows when a token is added to an already styled className", () => {
    const before = countAppearance("pages/X.tsx", '<i className="text-ink" />');
    const after = countAppearance("pages/X.tsx", '<i className="text-ink font-mono" />');
    expect(after).toBeGreaterThan(before);
  });

  it("ignores primitives, render layers and stories", () => {
    const src = '<span className="text-caption text-ink" />';
    expect(countAppearance("ui/Text.tsx", src)).toBe(0);
    expect(countAppearance("plot/TracePlot.tsx", src)).toBe(0);
    expect(countAppearance("pages/X.stories.tsx", src)).toBe(0);
  });

  it("counts opacity-slash and important forms", () => {
    expect(countAppearance("pages/X.tsx", '<i className="bg-ink/50 !font-mono hover:!text-ink/70" />')).toBe(3);
  });

  it("treats text alignment and wrapping as placement", () => {
    expect(countAppearance("pages/X.tsx", '<p className="text-center text-nowrap" />')).toBe(0);
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

  it("counts trailing and mid-line comments", () => {
    expect(countTsComments("const x = 1; // note")).toBe(1);
    expect(countTsComments("<a>{/* hidden */}</a>")).toBe(1);
  });

  it("does not count comment-like text inside strings", () => {
    expect(countTsComments('const url = "http://example.com"; const g = "/* no */";')).toBe(0);
  });

  it("exempts triple-slash reference directives", () => {
    expect(countTsComments('/// <reference types="vite/client" />')).toBe(0);
  });

  it("exempts directives with a reason and leading physics tags", () => {
    const src = [
      "// eslint-disable-next-line no-console -- surfaced to the user as a toast",
      "// @ts-expect-error upstream types",
      "// physics: q = 4π sin(θ)/λ",
    ].join("\n");
    expect(countTsComments(src)).toBe(0);
  });

  it("counts reason-less directives, @ts-ignore and a physics word mid-comment", () => {
    expect(countTsComments("// eslint-disable-next-line no-console")).toBe(1);
    expect(countTsComments("// @ts-ignore")).toBe(1);
    expect(countTsComments("// see physics: later")).toBe(1);
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

  it("counts indented docstrings", () => {
    const src = ["module M", '    """', "    Doc.", '    """', "    function g() end", "end"].join("\n");
    expect(countJuliaComments(src)).toBe(3);
  });

  it("counts trailing # comments and #= =# blocks", () => {
    expect(countJuliaComments("x = 1  # note")).toBe(1);
    expect(countJuliaComments(["#=", "block", "=#", "y = 2"].join("\n"))).toBe(3);
  });

  it("does not count # inside strings or chars, or SQL strings", () => {
    const src = [
      'msg = "item #3"',
      "c = '#'",
      "rows = DBInterface.execute(db,",
      '    """SELECT id',
      "       FROM t",
      '    """)',
      'q = """',
      "SELECT 1 -- # not a comment",
      '""")',
    ].join("\n");
    expect(countJuliaComments(src)).toBe(0);
  });

  it("exempts leading physics tags only", () => {
    expect(countJuliaComments("# physics: Q0f31 counts are log-normalised")).toBe(0);
    expect(countJuliaComments("x = 1  # physics: beam centre in pixels")).toBe(0);
    expect(countJuliaComments("# note on physics: later")).toBe(1);
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

describe("canUpdate", () => {
  it("allows a rename or split when each kind's total does not grow", () => {
    const baseline = { appearance: { "a.tsx": 3 }, comments: { "a.tsx": 10 } };
    const moved = { appearance: { "b.tsx": 3 }, comments: { "a.tsx": 4, "c.tsx": 6 } };
    expect(canUpdate(baseline, moved)).toBe(true);
  });

  it("refuses when a kind's total grows", () => {
    const baseline = { appearance: { "a.tsx": 3 }, comments: {} };
    expect(canUpdate(baseline, { appearance: { "a.tsx": 2, "b.tsx": 2 }, comments: {} })).toBe(false);
  });
});

describe("comment scope", () => {
  it("covers every Julia and TypeScript package, tests and tools included", () => {
    expect(JULIA_ROOTS).toEqual(
      expect.arrayContaining([
        "src",
        "test",
        "packages/HimalayaUI/src",
        "packages/HimalayaUI/test",
        "packages/HimalayaPeakModel",
        "packages/HimalayaDB",
        "tools",
      ]),
    );
    expect(TS_ROOTS).toEqual(
      expect.arrayContaining([
        "packages/HimalayaUI/frontend/src",
        "packages/HimalayaUI/frontend/test",
        "packages/HimalayaUI/frontend/e2e",
        "packages/HimalayaUI/frontend/scripts",
      ]),
    );
  });
});
