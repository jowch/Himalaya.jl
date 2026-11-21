import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, resolve } from "node:path";

const here = dirname(fileURLToPath(import.meta.url));
const styles = readFileSync(resolve(here, "../src/styles.css"), "utf8");
const fontsource = (w: number) => readFileSync(resolve(here, `../node_modules/@fontsource/source-sans-3/${w}.css`), "utf8");

interface Face { file: string; range: string; override: boolean }

function faces(css: string): Face[] {
  return [...css.matchAll(/@font-face\s*\{([^}]*)\}/g)].map(([, body]) => ({
    file: /files\/([\w-]+)\.woff2/.exec(body!)?.[1] ?? "",
    range: /unicode-range:\s*([^;]+);/.exec(body!)?.[1]?.trim() ?? "",
    override: /ascent-override/.test(body!) && /descent-override/.test(body!),
  }));
}

describe("Source Sans 3 metric override", () => {
  const overrides = faces(styles).filter((f) => f.override);
  for (const w of [400, 600]) {
    for (const subset of ["latin", "latin-ext"]) {
      it(`redeclares the ${w} ${subset} face with fontsource's exact unicode-range`, () => {
        const file = `source-sans-3-${subset}-${w}-normal`;
        const upstream = faces(fontsource(w)).find((f) => f.file === file);
        expect(upstream, `fontsource no longer ships ${file}`).toBeDefined();
        expect(overrides.find((f) => f.file === file)?.range).toBe(upstream!.range);
      });
    }
  }

  it("declares the overriding faces after the fontsource imports, so they win", () => {
    const lastImport = Math.max(styles.lastIndexOf('@fontsource/source-sans-3/400.css";'), styles.lastIndexOf('@fontsource/source-sans-3/600.css";'));
    expect(lastImport).toBeGreaterThan(-1);
    expect(styles.indexOf("ascent-override")).toBeGreaterThan(lastImport);
  });
});
