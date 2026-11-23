import { describe, it, expect } from "vitest";
import { cssBlock, stylesCss, THEME_BLOCKS } from "../helpers/cssBlocks";

const grow = (selector: string): number => Number(cssBlock(selector).get("--peak-manual-grow")?.match(/^([\d.]+)px$/)?.[1] ?? NaN);

describe("manual peak marker outline", () => {
  it("is drawn only in dark themes, where the light auto marker reads larger", () => {
    expect(grow(":root")).toBe(0);
    expect(grow(THEME_BLOCKS.system)).toBeGreaterThan(0);
    expect(grow(THEME_BLOCKS.dark)).toBe(grow(THEME_BLOCKS.system));
    expect(stylesCss).toMatch(/\.peak-manual\s*\{[^}]*stroke-width:\s*var\(--peak-manual-grow\)/);
  });
});
