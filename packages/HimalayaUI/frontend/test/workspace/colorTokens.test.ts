import { describe, it, expect } from "vitest";
import { readFileSync, readdirSync, statSync } from "node:fs";
import { join, resolve } from "node:path";

const SRC = resolve(__dirname, "../../src");

function files(dir: string): string[] {
  return readdirSync(dir).flatMap((name) => {
    const p = join(dir, name);
    return statSync(p).isDirectory() ? files(p) : /\.(ts|tsx)$/.test(name) ? [p] : [];
  });
}

function outsideTheme(css: string): string {
  const start = css.indexOf("@theme");
  if (start < 0) return css;
  let depth = 0;
  for (let i = css.indexOf("{", start); i < css.length; i++) {
    if (css[i] === "{") depth++;
    if (css[i] === "}" && --depth === 0) return css.slice(0, start) + css.slice(i + 1);
  }
  return css;
}

describe("colour tokens the workspace reads from script", () => {
  it("are emitted by Tailwind, so the browser can resolve them", () => {
    const css = readFileSync(join(SRC, "styles.css"), "utf8");
    const plain = outsideTheme(css);
    const all = files(SRC).map((f) => readFileSync(f, "utf8")).join("\n");
    const used = new Set<string>();
    for (const f of files(join(SRC, "workspace"))) {
      for (const m of readFileSync(f, "utf8").matchAll(/var\(--color-([a-z0-9-]+)\)/g)) used.add(m[1]!);
    }
    for (const m of readFileSync(join(SRC, "ui/WorkspacePlot.tsx"), "utf8").matchAll(/--color-([a-z]+)-\$\{n\}/g)) {
      for (let n = 1; n <= 5; n++) used.add(`${m[1]}-${n}`);
    }
    const missing = [...used].filter((name) => {
      if (new RegExp(`--color-${name}\\s*:`).test(plain)) return false;
      if (!new RegExp(`--color-${name}\\s*:`).test(css)) return true;
      return !new RegExp(`\\b(bg|text|fill|stroke|border|outline|ring|shadow|decoration|from|to|via)-${name}\\b`).test(all);
    });
    expect(missing).toEqual([]);
  });
});
