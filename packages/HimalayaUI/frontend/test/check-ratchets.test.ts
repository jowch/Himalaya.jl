import { execFileSync } from "node:child_process";
import { describe, it, expect } from "vitest";
import {
  countAppearance,
  countTsComments,
  countJuliaComments,
  compareToBaseline,
  canUpdate,
  countCssComments,
  trackedSources,
  measure,
  EXCLUDED,
  DIRECTIVE_MAX,
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

  it("counts one-line, @doc and raw docstrings", () => {
    expect(countJuliaComments(['"Doc."', "f(x) = 1"].join("\n"))).toBe(1);
    expect(countJuliaComments(['"Doc."', "function f(x)", "end"].join("\n"))).toBe(1);
    expect(countJuliaComments(['"Doc."', "const LIMIT = 3"].join("\n"))).toBe(1);
    expect(countJuliaComments(['"Doc."', "first_peak_q(::Type{Lamellar}, a) = 2π / a"].join("\n"))).toBe(1);
    expect(countJuliaComments(['@doc """', "Doc.", '""" f'].join("\n"))).toBe(3);
    expect(countJuliaComments('@doc "Doc." f')).toBe(1);
    expect(countJuliaComments(['raw"""', "Doc \\alpha.", '"""', "f(x) = 1"].join("\n"))).toBe(3);
  });

  it("counts one-line docstrings on fields, operator methods and callable structs", () => {
    expect(countJuliaComments(["struct S", '    "doc"', "    x::Int", "end"].join("\n"))).toBe(1);
    expect(countJuliaComments(["Base.@kwdef struct C", '    "doc"', "    x::Float64 = 1.0", "end"].join("\n"))).toBe(1);
    expect(countJuliaComments(['"doc"', "Base.:+(a::S, b::S) = a"].join("\n"))).toBe(1);
    expect(countJuliaComments(['"doc"', "(s::S)(x) = 1"].join("\n"))).toBe(1);
  });

  it("does not count a string before a call or a comparison", () => {
    expect(countJuliaComments(["begin", '    "x"', "    push!(v, k = 1)", "end"].join("\n"))).toBe(0);
    expect(countJuliaComments(["begin", '    "x"', "    f(a) == b", "end"].join("\n"))).toBe(0);
  });

  it("does not count a string returned at the end of a block", () => {
    expect(countJuliaComments(["function f()", '    "done"', "end"].join("\n"))).toBe(0);
    expect(countJuliaComments(['x = f(', '    "arg"', ")"].join("\n"))).toBe(0);
  });

  it("does not count # inside interpolated quotes or command literals", () => {
    expect(countJuliaComments('s = "a $(f("#")) b"')).toBe(0);
    expect(countJuliaComments("run(`echo # x`)")).toBe(0);
    expect(countJuliaComments('m = r"#\\d+"')).toBe(0);
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
  const gone = () => false;

  it("accepts counts that only fall", () => {
    const baseline = { appearance: { "a.tsx": 3 }, comments: { "a.tsx": 10 } };
    expect(canUpdate(baseline, { appearance: { "a.tsx": 2 }, comments: {} }, gone)).toBe(true);
  });

  it("refuses moving allowance from one file to another", () => {
    const baseline = { appearance: {}, comments: { "a.jl": 10, "b.jl": 2 } };
    expect(canUpdate(baseline, { appearance: {}, comments: { "a.jl": 4, "b.jl": 8 } }, gone)).toBe(false);
  });

  it("refuses a new file with any count", () => {
    const baseline = { appearance: {}, comments: { "a.jl": 10 } };
    expect(canUpdate(baseline, { appearance: {}, comments: { "a.jl": 9, "c.jl": 1 } }, gone)).toBe(false);
  });

  it("accepts a pure rename to a file with the same count", () => {
    const baseline = { appearance: { "a.tsx": 3 }, comments: { "a.tsx": 10 } };
    expect(canUpdate(baseline, { appearance: { "b.tsx": 3 }, comments: { "b.tsx": 10 } }, gone)).toBe(true);
  });

  it("refuses a split, and a rename whose old file still exists", () => {
    const baseline = { appearance: {}, comments: { "a.tsx": 10 } };
    expect(canUpdate(baseline, { appearance: {}, comments: { "b.tsx": 4, "c.tsx": 6 } }, gone)).toBe(false);
    expect(canUpdate(baseline, { appearance: {}, comments: { "b.tsx": 10 } }, () => true)).toBe(false);
  });
});

describe("countCssComments", () => {
  it("counts block comment lines and ignores strings", () => {
    expect(countCssComments(["/* one", "   two */", "a { color: red; }"].join("\n"))).toBe(2);
    expect(countCssComments('a::before { content: "/* no */"; }')).toBe(0);
    expect(countCssComments("/* physics: q = 4π sin(θ)/λ */")).toBe(0);
  });
});

describe("directive length", () => {
  it("counts a directive whose reason runs past the cap or onto more lines", () => {
    const long = `// @ts-expect-error ${"x".repeat(DIRECTIVE_MAX)}`;
    expect(countTsComments(long)).toBe(1);
    expect(countTsComments("/* @ts-expect-error a\n   b */")).toBe(2);
  });
});

describe("comment scope", () => {
  const SOURCE = /\.(?:jl|[cm]?[jt]sx?|css)$/;
  const allTracked = execFileSync("git", ["ls-files", "--full-name", ":/"], { encoding: "utf8" }).split("\n");

  it("measures every tracked Julia, TypeScript, JavaScript and CSS file unless excluded", () => {
    const measured = new Set(trackedSources());
    const missing = allTracked.filter((f) => SOURCE.test(f) && !measured.has(f) && !EXCLUDED.includes(f));
    expect(missing).toEqual([]);
  });

  it("names only tracked files in the exclude list", () => {
    expect(EXCLUDED.filter((f) => !allTracked.includes(f))).toEqual([]);
  });

  it("counts comments in stories, config files and every package", () => {
    const keys = Object.keys(measure().comments ?? {});
    for (const prefix of ["examples/", "packages/HimalayaUI/bench/", "packages/HimalayaUI/frontend/vite.config.ts", "packages/HimalayaUI/frontend/.storybook/", "src/", "tools/"]) {
      expect(keys.some((k) => k.startsWith(prefix))).toBe(true);
    }
    expect(keys.some((k) => k.includes(".stories."))).toBe(true);
  });
});
