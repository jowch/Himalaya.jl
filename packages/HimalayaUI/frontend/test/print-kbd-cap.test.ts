import { describe, it, expect } from "vitest";
import { readdirSync, readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join, relative, resolve } from "node:path";

/*
 * Keyboard-hint law: inline key hints use the KbKey cap, not dimmed glyphs.
 * ---------------------------------------------------------------------------
 * Bare `font-mono text-xs opacity-60` glyphs (e.g. Prev`[`, Drop`X`) composite
 * their inherited ink down to ~2.7:1 — sub-AA — AND diverge from the boxed
 * KbKey caps the KbdLegend uses (CC-KBDCAP). Every inline shortcut hint must
 * render through the <KbKey> primitive (bg-plate / text-ink-soft, 6.8:1).
 *
 * The scan bans the `opacity-60 ml-1` inline-key-hint idiom anywhere under
 * src/print. Mirrors test/print-copy-no-emdash.test.ts (TEST-ONLY, reads the
 * real sources at test time). ui/ is excluded: the primitives own appearance.
 */

const __dirname = dirname(fileURLToPath(import.meta.url));
const PRINT_ROOT = resolve(__dirname, "../src/print");
const UI_DIR = resolve(PRINT_ROOT, "ui");

function walk(dir: string, out: string[] = []): string[] {
  for (const entry of readdirSync(dir, { withFileTypes: true })) {
    const p = join(dir, entry.name);
    if (p.startsWith(UI_DIR)) continue;
    if (entry.isDirectory()) walk(p, out);
    else if (/\.tsx?$/.test(entry.name) && !/\.stories\.tsx?$/.test(entry.name)) out.push(p);
  }
  return out;
}

describe("keyboard-hint law — src/print inline key hints", () => {
  it("renders inline key hints via KbKey, not dimmed bare glyphs (CC-KBDCAP)", () => {
    const offenders: string[] = [];
    for (const file of walk(PRINT_ROOT)) {
      readFileSync(file, "utf8")
        .split("\n")
        .forEach((line, i) => {
          if (/opacity-60\s+ml-1\b/.test(line)) {
            offenders.push(`${relative(PRINT_ROOT, file)}:${i + 1}: ${line.trim()}`);
          }
        });
    }
    expect(offenders).toEqual([]);
  });
});
