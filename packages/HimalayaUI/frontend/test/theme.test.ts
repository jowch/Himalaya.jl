import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { createRequire } from "node:module";
import { resolve } from "node:path";
import { compile } from "tailwindcss";

const require = createRequire(import.meta.url);
const STYLES_PATH = resolve(__dirname, "../src/styles.css");
const SRC = readFileSync(STYLES_PATH, "utf-8");

// Extract the contents of the top-level `@theme { ... }` block from source.
// (Mirrors test/phases.test.ts: assert against source text, not a runtime DOM.)
function themeBlock(src: string): string {
  const open = src.indexOf("@theme");
  if (open < 0) throw new Error("no @theme block");
  const brace = src.indexOf("{", open);
  let depth = 0;
  for (let i = brace; i < src.length; i++) {
    if (src[i] === "{") depth++;
    else if (src[i] === "}") {
      depth--;
      if (depth === 0) return src.slice(brace + 1, i);
    }
  }
  throw new Error("unterminated @theme block");
}

// Read a `--name: value;` declaration's value from a CSS text region.
function declValue(css: string, name: string): string | null {
  const m = new RegExp(`${name.replace(/[-]/g, "\\-")}\\s*:\\s*([^;]+);`).exec(css);
  return m ? m[1]!.trim() : null;
}

// Compile the real stylesheet and return the resolved value of the var chain a
// utility's declared property points at (e.g. `border-radius: var(--radius-md)`
// resolved against the emitted `--radius-md: 5px`).
async function resolvedUtilityValue(utility: string, property: string): Promise<string> {
  const compiler = await compile(SRC, {
    base: resolve(__dirname, ".."),
    loadStylesheet: async (id: string, base: string) => {
      if (id === "tailwindcss") {
        const p = require.resolve("tailwindcss/index.css");
        return { base, content: readFileSync(p, "utf-8"), path: p };
      }
      // fontsource + any other import: irrelevant to token resolution
      return { base, content: "", path: id };
    },
  });
  const out = compiler.build([utility]);
  const ruleRe = new RegExp(`\\.${utility.replace(/[\\^$.*+?()[\]{}|]/g, "\\$&")}\\s*\\{([^}]*)\\}`);
  const rule = ruleRe.exec(out);
  if (!rule) throw new Error(`utility .${utility} not generated`);
  const propRe = new RegExp(`${property}\\s*:\\s*([^;]+);`);
  const prop = propRe.exec(rule[1]!);
  if (!prop) throw new Error(`${property} not set by .${utility}`);
  let value = prop[1]!.trim();
  // Follow a single var(--x) hop to its emitted :root declaration.
  const varRef = /^var\((--[\w-]+)\)$/.exec(value);
  if (varRef) {
    const resolved = declValue(out, varRef[1]!);
    if (!resolved) throw new Error(`${varRef[1]} not emitted in compiled output`);
    value = resolved;
  }
  return value;
}

describe("styles.css @theme tokens (Phase 0 token foundation)", () => {
  const theme = themeBlock(SRC);

  it("defines --color-scrim", () => {
    expect(declValue(theme, "--color-scrim")).toBe("oklch(0.05 0 0 / 0.65)");
  });

  it("sources --color-print-accent from --color-accent (single source, no duplicated literal)", () => {
    expect(declValue(theme, "--color-print-accent")).toBe("var(--color-accent)");
  });

  it("keeps the singular --radius (= 5px) so bare `rounded` snaps to 5px, not stock 4px", () => {
    // Tailwind v4's bare `.rounded` utility emits `border-radius: var(--radius)`.
    // Deleting --radius silently regresses every bare `rounded` to the stock
    // 0.25rem fallback — so it must stay defined and on the 5px control scale.
    expect(declValue(theme, "--radius")).toBe("5px");
  });

  it("defines the namespaced radius scale", () => {
    expect(declValue(theme, "--radius-sm")).toBe("5px");
    expect(declValue(theme, "--radius-md")).toBe("5px");
    expect(declValue(theme, "--radius-full")).toBe("9999px");
  });

  it("resolves rounded-md utility to 5px (not stock Tailwind 6px)", async () => {
    expect(await resolvedUtilityValue("rounded-md", "border-radius")).toBe("5px");
  });

  it("resolves rounded-sm utility to 5px (not stock Tailwind 4px)", async () => {
    expect(await resolvedUtilityValue("rounded-sm", "border-radius")).toBe("5px");
  });

  it("resolves bare `rounded` utility to 5px (regression guard: not stock Tailwind 4px)", async () => {
    expect(await resolvedUtilityValue("rounded", "border-radius")).toBe("5px");
  });

  it("no longer hardcodes a 12px corner on .card", () => {
    const cardRule = /\.card\s*\{([^}]*)\}/.exec(SRC);
    expect(cardRule).not.toBeNull();
    expect(/border-radius\s*:\s*12px/.test(cardRule![1]!)).toBe(false);
  });
});
