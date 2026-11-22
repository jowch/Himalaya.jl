import { describe, it, expect } from "vitest";
import { readdirSync, readFileSync, statSync } from "node:fs";
import { join, relative } from "node:path";

const SRC = join(__dirname, "../../src");
const ROOTS = [join(SRC, "workspace"), join(SRC, "ui")];
const PAGE_LAYERS = /from\s+["']\.\.\/(pages|components|shell)\//;

function files(dir: string): string[] {
  return readdirSync(dir).flatMap((n) => {
    const p = join(dir, n);
    return statSync(p).isDirectory() ? files(p) : /\.(ts|tsx)$/.test(n) ? [p] : [];
  });
}

describe("workspace and ui import boundary", () => {
  it("import nothing from the pages, components or shell layers", () => {
    const offenders = ROOTS.flatMap(files).flatMap((f) =>
      readFileSync(f, "utf8").split("\n").filter((l) => PAGE_LAYERS.test(l)).map((l) => `${relative(SRC, f)}: ${l.trim()}`),
    );
    expect(offenders).toEqual([]);
  });

  it("would catch a page-layer import", () => {
    expect(PAGE_LAYERS.test('import { X } from "../shell/X";')).toBe(true);
  });
});
