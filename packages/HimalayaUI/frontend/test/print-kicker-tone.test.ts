import { describe, it, expect } from "vitest";
import { readdirSync, readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join, relative, resolve } from "node:path";

/*
 * Kicker tone law: no `tone="faint"` on content labels (CC-KICKER-FAINT).
 * ---------------------------------------------------------------------------
 * `ink-faint` measures 3.16:1 on paper / 3.29:1 on plate / 2.92:1 on the sunk
 * rail (see test/contrast-tokens.test.ts + styles.css:343-351). Kickers are
 * 700-weight uppercase at text-sm (~13px) — NOT WCAG "large text" — so the
 * 4.5:1 AA-normal line applies and `faint` fails it on every surface. The
 * documented rule reserves `faint` for purely decorative labels on plain
 * paper; every section eyebrow / panel header / `h2` that LABELS content must
 * use `tone="soft"` (ink-soft, AA at all sizes).
 *
 * No call site currently has a legitimate decorative-on-paper use, so the
 * floor is zero. If a genuinely decorative kicker ever appears, add its
 * `path:line` to ALLOWLIST with a one-line justification.
 *
 * Mirrors the file-reading guard pattern of test/print-copy-no-emdash.test.ts:
 * TEST-ONLY, reads the real sources at test time. The Kicker primitive's own
 * `tone = "faint"` DEFAULT (a destructured param, written `tone = "faint"`
 * with spaces) is not a JSX prop and never matches the `tone="faint"` scan;
 * src/print/ui is excluded regardless.
 */

const __dirname = dirname(fileURLToPath(import.meta.url));
const PRINT_ROOT = resolve(__dirname, "../src/print");
const UI_DIR = resolve(PRINT_ROOT, "ui");

/** path:line entries that are genuinely decorative on plain paper. Empty today. */
const ALLOWLIST = new Set<string>([]);

function walk(dir: string, out: string[] = []): string[] {
  for (const entry of readdirSync(dir, { withFileTypes: true })) {
    const p = join(dir, entry.name);
    if (p.startsWith(UI_DIR)) continue; // primitives own the appearance, incl. the default tone
    if (entry.isDirectory()) walk(p, out);
    else if (/\.tsx?$/.test(entry.name) && !/\.stories\.tsx?$/.test(entry.name)) out.push(p);
  }
  return out;
}

describe("kicker tone law — src/print content labels", () => {
  it('uses no tone="faint" outside the decorative allowlist (CC-KICKER-FAINT)', () => {
    const offenders: string[] = [];
    for (const file of walk(PRINT_ROOT)) {
      readFileSync(file, "utf8")
        .split("\n")
        .forEach((line, i) => {
          // Both forms: the Kicker primitive's `tone="faint"` AND a raw
          // `text-kicker-faint` utility hand-applied to a content label.
          if (line.includes('tone="faint"') || line.includes("text-kicker-faint")) {
            const ref = `${relative(PRINT_ROOT, file)}:${i + 1}`;
            if (!ALLOWLIST.has(ref)) offenders.push(`${ref}: ${line.trim()}`);
          }
        });
    }
    expect(offenders).toEqual([]);
  });
});
