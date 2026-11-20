import { describe, it, expect } from "vitest";
import { readdirSync, readFileSync, statSync } from "node:fs";
import { join, resolve } from "node:path";
import { effectiveTheme } from "../../src/workspace/theme";

const css = readFileSync(resolve(__dirname, "../../src/styles.css"), "utf8");

function block(selector: string): Map<string, string> {
  const start = css.indexOf(selector + " {");
  expect(start).toBeGreaterThan(-1);
  const body = css.slice(css.indexOf("{", start) + 1, css.indexOf("}", start));
  return new Map(
    body
      .split(";")
      .map((d) => d.trim())
      .filter((d) => d.startsWith("--"))
      .map((d) => [d.slice(0, d.indexOf(":")).trim(), d.slice(d.indexOf(":") + 1).trim()] as [string, string]),
  );
}

describe("workspace theme", () => {
  it("follows the system setting unless the user picked one", () => {
    expect(effectiveTheme("system", true)).toBe("dark");
    expect(effectiveTheme("system", false)).toBe("light");
    expect(effectiveTheme("light", true)).toBe("light");
    expect(effectiveTheme("dark", false)).toBe("dark");
  });

  it("states the same dark tokens for the system setting and the explicit choice", () => {
    expect(block('[data-workspace-theme="system"]')).toEqual(block('[data-workspace-theme="dark"]'));
  });

  it("restates every token that mixes another token, since those resolve where they are declared", () => {
    const dark = block('[data-workspace-theme="dark"]');
    const theme = css.slice(css.indexOf("@theme {"), css.indexOf("\n}\n", css.indexOf("@theme {")));
    const derivedColours = Array.from(theme.matchAll(/^\s*(--color-[\w-]+):[^;]*var\(--color-/gm)).map((m) => m[1]!);
    const files = (dir: string): string[] =>
      readdirSync(dir).flatMap((f) => (statSync(join(dir, f)).isDirectory() ? files(join(dir, f)) : [join(dir, f)]));
    const source = files(resolve(__dirname, "../../src/workspace")).map((f) => readFileSync(f, "utf8")).join("\n");
    const used = derivedColours.filter((name) => {
      const short = name.replace("--color-", "");
      return source.includes(name) || new RegExp(`-${short}\\b`).test(source);
    });
    expect(used).toContain("--color-row-cursor");
    for (const name of used) expect([name, dark.has(name)]).toEqual([name, true]);
  });
});
