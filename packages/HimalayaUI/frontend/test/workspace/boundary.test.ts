import { describe, it, expect } from "vitest";
import { readdirSync, readFileSync, statSync } from "node:fs";
import { join, relative } from "node:path";

const ROOT = join(__dirname, "../../src/workspace");
const OLD = /from\s+["'][./]*\/print\/(pages|components|shell|interaction)\//;

function files(dir: string): string[] {
  return readdirSync(dir).flatMap((n) => {
    const p = join(dir, n);
    return statSync(p).isDirectory() ? files(p) : /\.(ts|tsx)$/.test(n) ? [p] : [];
  });
}

describe("workspace import boundary", () => {
  it("imports nothing from the old pages, components, shell or interaction layers", () => {
    const offenders = files(ROOT).flatMap((f) =>
      readFileSync(f, "utf8").split("\n").filter((l) => OLD.test(l)).map((l) => `${relative(ROOT, f)}: ${l.trim()}`),
    );
    expect(offenders).toEqual([]);
  });
});
