import { describe, it, expect } from "vitest";
import { readdirSync, readFileSync, statSync } from "node:fs";
import { dirname, join, relative, resolve } from "node:path";

const SRC = join(__dirname, "../../src");
const ROOTS = [join(SRC, "workspace"), join(SRC, "ui")];
const PAGE_LAYERS = ["pages", "components", "shell"];
const RELATIVE_IMPORT = /from\s+["'](\.{1,2}\/[^"']+)["']/g;

function files(dir: string): string[] {
  return readdirSync(dir).flatMap((n) => {
    const p = join(dir, n);
    return statSync(p).isDirectory() ? files(p) : /\.(ts|tsx)$/.test(n) ? [p] : [];
  });
}

function pageLayerImports(file: string, content: string): string[] {
  return [...content.matchAll(RELATIVE_IMPORT)]
    .map((m) => m[1]!)
    .filter((spec) => PAGE_LAYERS.includes(relative(SRC, resolve(dirname(file), spec)).split(/[\\/]/)[0]!));
}

describe("workspace and ui import boundary", () => {
  it("import nothing from the pages, components or shell layers", () => {
    const offenders = ROOTS.flatMap(files).flatMap((f) =>
      pageLayerImports(f, readFileSync(f, "utf8")).map((spec) => `${relative(SRC, f)}: ${spec}`),
    );
    expect(offenders).toEqual([]);
  });

  it("catches a page-layer import at any depth", () => {
    const deep = join(SRC, "workspace", "a", "b", "X.tsx");
    expect(pageLayerImports(deep, 'import { X } from "../../../shell/X";')).toEqual(["../../../shell/X"]);
    expect(pageLayerImports(join(SRC, "ui", "X.tsx"), 'import { X } from "../components/X";')).toEqual(["../components/X"]);
    expect(pageLayerImports(join(SRC, "ui", "X.tsx"), 'import { X } from "../lib/X";')).toEqual([]);
  });
});
