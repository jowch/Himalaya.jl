import { describe, it, expect } from "vitest";
import { readdirSync, readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join, relative, resolve } from "node:path";

/*
 * House copy law: NO EM DASHES in user-facing copy.
 * ---------------------------------------------------------------------------
 * This pin comment-strips every non-story source under src/print and asserts
 * no "—" survives outside the one sanctioned use: the standalone "—" string
 * literal used as a missing-value placeholder (e.g. `?? "—"`, `value: "—"`).
 *
 * Mirrors the file-reading guard pattern of test/contrast-tokens.test.ts:
 * TEST-ONLY, reads the real sources at test time.
 *
 * Deliberate scope limits:
 *  - Code COMMENTS may use em dashes (stripped before scanning).
 *  - *.stories.tsx are dev-facing narration, excluded from the zero criterion
 *    (stories that mirror page copy are kept in sync by hand).
 *  - "--" is NOT checked here: CSS custom-property names (`--color-*`) appear
 *    legitimately in string literals, so a "--" scan cannot discriminate.
 *  - Comment stripping is regex-naive (a "//" inside a string truncates that
 *    line's scan). That can only under-report, never false-positive.
 */

const __dirname = dirname(fileURLToPath(import.meta.url));
const PRINT_ROOT = resolve(__dirname, "../src/print");

function walk(dir: string, out: string[] = []): string[] {
  for (const entry of readdirSync(dir, { withFileTypes: true })) {
    const p = join(dir, entry.name);
    if (entry.isDirectory()) walk(p, out);
    else if (/\.tsx?$/.test(entry.name) && !/\.stories\.tsx?$/.test(entry.name)) out.push(p);
  }
  return out;
}

/** Blank out block + line comments, preserving line numbers. */
function stripComments(src: string): string {
  return src
    .replace(/\/\*[\s\S]*?\*\//g, (m) => m.replace(/[^\n]/g, " "))
    .replace(/\/\/[^\n]*/g, "");
}

/** The sanctioned placeholder: a string literal that is EXACTLY one em dash. */
const PLACEHOLDER_LITERAL = /(["'`])—\1/g;

describe("house copy law — src/print prose", () => {
  it("contains no em dash outside comments and standalone placeholder literals", () => {
    const offenders: string[] = [];
    for (const file of walk(PRINT_ROOT)) {
      const stripped = stripComments(readFileSync(file, "utf8"));
      stripped.split("\n").forEach((line, i) => {
        if (line.replace(PLACEHOLDER_LITERAL, "").includes("—")) {
          offenders.push(`${relative(PRINT_ROOT, file)}:${i + 1}: ${line.trim()}`);
        }
      });
    }
    expect(offenders).toEqual([]);
  });
});
