import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const css = readFileSync(resolve(__dirname, "../../src/styles.css"), "utf8");

function block(selector: string): string {
  const at = css.indexOf(selector);
  expect(at).toBeGreaterThanOrEqual(0);
  return css.slice(at, css.indexOf("}", at));
}
const grow = (b: string): number => Number(b.match(/--peak-manual-grow:\s*([\d.]+)px/)?.[1] ?? NaN);

describe("manual peak marker outline", () => {
  it("is drawn only in dark themes, where the light auto marker reads larger", () => {
    expect(grow(block(":root {"))).toBe(0);
    expect(grow(block('[data-theme="system"] {'))).toBeGreaterThan(0);
    expect(grow(block('[data-theme="dark"] {'))).toBe(grow(block('[data-theme="system"] {')));
    expect(css).toMatch(/\.peak-manual\s*\{[^}]*stroke-width:\s*var\(--peak-manual-grow\)/);
  });
});
