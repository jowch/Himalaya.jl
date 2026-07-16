import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, resolve } from "node:path";

/*
 * F-STATE focus-visible pin — guards the systematized keyboard-focus state.
 * ---------------------------------------------------------------------------
 * DESIGN.md (State taxonomy → Focus-visible): every bespoke <button>/<a>
 * inherits the terracotta accent ring from ONE base-layer rule in styles.css,
 * instead of ~15 per-site class edits (or, worse, the UA-default ring).
 * Primitives that refine the ring (ghost 1px, etc.) win via the utilities
 * layer; deliberate suppressions win via @layer components or focus:outline-*
 * utilities — Tailwind v4 layer order is base < components < utilities.
 *
 * Like contrast-tokens.test.ts, this reads the REAL stylesheet text at test
 * time. It pins the stylesheet itself (the allowed pattern), not Tailwind
 * class strings on components.
 */

const __dirname = dirname(fileURLToPath(import.meta.url));
const STYLES_CSS = resolve(__dirname, "../src/styles.css");
const css = readFileSync(STYLES_CSS, "utf8");

// Extract the full `@layer base { … }` block by brace counting (it nests
// keyframes/rules, so a lazy regex would stop at the first `}`).
function extractLayerBase(source: string): string {
  const start = source.search(/@layer\s+base\s*\{/);
  if (start === -1) return "";
  const open = source.indexOf("{", start);
  let depth = 0;
  for (let i = open; i < source.length; i++) {
    if (source[i] === "{") depth++;
    else if (source[i] === "}") {
      depth--;
      if (depth === 0) return source.slice(open + 1, i);
    }
  }
  return "";
}

const layerBase = extractLayerBase(css);

describe("focus-visible base rule — F-STATE (styles.css)", () => {
  it("styles.css has an @layer base block", () => {
    expect(layerBase.length).toBeGreaterThan(0);
  });

  it("@layer base contains a button:focus-visible + a:focus-visible rule", () => {
    expect(layerBase).toMatch(
      /button:focus-visible\s*,\s*a:focus-visible\s*\{/,
    );
  });

  it("the rule sets the 2px accent outline with 2px offset", () => {
    const m = layerBase.match(
      /button:focus-visible\s*,\s*a:focus-visible\s*\{([^}]*)\}/,
    );
    expect(m, "button/a focus-visible rule body present").toBeTruthy();
    const body = m![1]!;
    expect(body).toMatch(/outline:\s*2px\s+solid\s+var\(--color-accent\)\s*;/);
    expect(body).toMatch(/outline-offset:\s*2px\s*;/);
  });
});
