import { readFileSync } from "node:fs";
import { resolve } from "node:path";

export const stylesCss = readFileSync(resolve(__dirname, "../../src/styles.css"), "utf8").replace(/\/\*[\s\S]*?\*\//g, "");

export function cssBlock(selector: string, css: string = stylesCss): Map<string, string> {
  const start = css.indexOf(selector + " {");
  if (start < 0) throw new Error(`${selector} not found in styles.css`);
  const body = css.slice(css.indexOf("{", start) + 1, css.indexOf("}", start));
  return new Map(
    body
      .split(";")
      .map((d) => d.trim())
      .filter((d) => d.startsWith("--"))
      .map((d) => [d.slice(0, d.indexOf(":")).trim(), d.slice(d.indexOf(":") + 1).trim()] as [string, string]),
  );
}

export const THEME_BLOCKS = {
  light: "@theme",
  dark: '[data-theme="dark"]',
  system: '[data-theme="system"]',
} as const;

export function themeToken(theme: keyof typeof THEME_BLOCKS, name: string): string {
  const value = cssBlock(THEME_BLOCKS[theme]).get(`--color-${name}`);
  if (!value) throw new Error(`--color-${name} not found in the ${theme} block`);
  return value;
}
