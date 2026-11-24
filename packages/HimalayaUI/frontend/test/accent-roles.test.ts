import { describe, it, expect } from "vitest";
import { contrastRatio } from "./helpers/contrast";
import { themeToken } from "./helpers/cssBlocks";

const AA_TEXT = 4.5;
const AA_NON_TEXT = 3;

describe("accent roles", () => {
  it("light fills use the accent itself with paper text", () => {
    expect(themeToken("light", "accent-fill")).toBe("var(--color-accent)");
    expect(themeToken("light", "on-accent")).toBe("var(--color-paper)");
    expect(contrastRatio(themeToken("light", "accent"), themeToken("light", "paper"))).toBeGreaterThanOrEqual(AA_TEXT);
  });

  it("dark fills hold white text, which the brighter dark accent could not", () => {
    const fill = themeToken("dark", "accent-fill");
    const onFill = themeToken("dark", "on-accent");
    expect(contrastRatio(onFill, fill)).toBeGreaterThanOrEqual(AA_TEXT);
    expect(contrastRatio(onFill, themeToken("dark", "accent"))).toBeLessThan(AA_TEXT);
  });

  it("dark accent text reads on the dark page and plate", () => {
    const accent = themeToken("dark", "accent");
    expect(contrastRatio(accent, themeToken("dark", "paper"))).toBeGreaterThanOrEqual(AA_TEXT);
    expect(contrastRatio(accent, themeToken("dark", "plate"))).toBeGreaterThanOrEqual(AA_TEXT);
  });

  it("a dark fill still stands off the page as a shape", () => {
    expect(contrastRatio(themeToken("dark", "accent-fill"), themeToken("dark", "paper"))).toBeGreaterThanOrEqual(AA_NON_TEXT);
  });

  it("the selected row is a neutral lift that keeps ink readable", () => {
    for (const theme of ["light", "dark"] as const) {
      const row = themeToken(theme, "row-selected");
      expect(row).toMatch(/oklch\(\s*[\d.]+\s+0\.0[01]\d/);
      expect(contrastRatio(themeToken(theme, "ink"), row)).toBeGreaterThanOrEqual(7);
    }
  });
});
