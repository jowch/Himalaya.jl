#!/usr/bin/env node
// Zero-dependency design-token guard. Globs src/**/*.{ts,tsx}, excludes src/components/ui/**,
// flags banned appearance utilities / raw color literals against a content-hash baseline.
// Pure functions are exported for unit testing; the CLI runs only when invoked directly.
import { readdirSync, readFileSync, writeFileSync, existsSync } from "node:fs";
import { join, relative, sep } from "node:path";
import { fileURLToPath } from "node:url";
import { createHash } from "node:crypto";

const HERE = fileURLToPath(new URL(".", import.meta.url));
const SRC_DIR = join(HERE, "..", "src");
const BASELINE_PATH = join(HERE, "design-baseline.json");

// --- Ban rules ---------------------------------------------------------------
// Each rule: { id, test(line) -> matched-substring|null }. Rule #3 also consults `relPath`
// (POSIX-normalized, relative to src/) against its allowlist.
//
// Rule #3 allowlist: arbitrary-value color UTILITIES that legitimately need a runtime color.
// (spec §4: rules #1/#2/#4/#5 take NO allowlist.)
const RULE3_ALLOWLIST = new Set([
  "phases.ts",
  "lib/comparison/coloring.ts",
  "components/MemberHeatmapLayer.tsx",
  "components/DetectorImage.tsx",
  "components/FocusDetectorPanel.tsx",
  "main.tsx",
]);
// figure-export/** is allowlisted by prefix (the whole export palette dir).
const RULE3_ALLOW_PREFIXES = ["lib/figure-export/"];

function rule3Allowed(relPath) {
  if (RULE3_ALLOWLIST.has(relPath)) return true;
  return RULE3_ALLOW_PREFIXES.some((p) => relPath.startsWith(p));
}

const COLOR_LITERAL = "(?:oklch|oklab|hsla?|rgba?)\\(|#[0-9a-fA-F]{3,8}\\b";

const RULES = [
  {
    // #1 inline arbitrary type size or color: text-[10px], text-[var(--color-…)]
    id: "no-arbitrary-text",
    test: (line) => {
      const m = /\btext-\[/.exec(line);
      return m ? line.slice(m.index, m.index + 24) : null;
    },
  },
  {
    // #2 inline arbitrary radius: rounded-[3px], rounded-[1.5px], rounded-[7px], …
    id: "no-arbitrary-radius",
    test: (line) => {
      const m = /\brounded-\[/.exec(line);
      return m ? line.slice(m.index, m.index + 20) : null;
    },
  },
  {
    // #3 raw appearance color inside a color UTILITY (not shadow-[…]).
    // shadow-[…] is excluded — 4 legit Plate-Lift rgba shadows live in tsx.
    id: "no-raw-color-utility",
    allowlisted: rule3Allowed,
    test: (line) => {
      const re = new RegExp(
        "\\b(?:bg|text|border|fill|stroke|ring|outline|from|via|to|decoration)-\\[(?:[^\\]]*?)(?:" +
          COLOR_LITERAL +
          ")",
      );
      const m = re.exec(line);
      return m ? line.slice(m.index, m.index + 40) : null;
    },
  },
  {
    // #4 side-stripe: border-l-4 / border-r-4 / border-{l,r}-[…] / border-{l,r}-2..9
    id: "no-side-stripe",
    test: (line) => {
      const m = /\bborder-(?:l|r)-(?:4|[2-9]|\[)/.exec(line);
      return m ? line.slice(m.index, m.index + 16) : null;
    },
  },
  {
    // #5 raw color literal ANYWHERE on the line (covers multi-line style={{…}}).
    // Strip every var(--color-*) token FIRST, then flag a remaining raw literal.
    // This passes legit computed-color sites (their non-var color is a ${…} expression
    // with no raw literal) and flags the raw oklch(0.05…) scrim hand-inlines.
    id: "no-raw-color-literal",
    test: (line) => {
      const stripped = line.replace(/var\(--color-[^)]*\)/g, "");
      const re = new RegExp("(?:" + COLOR_LITERAL + ")");
      const m = re.exec(stripped);
      return m ? stripped.slice(m.index, m.index + 32).trim() : null;
    },
  },
];

// --- Violation discovery -----------------------------------------------------
function listSourceFiles(dir) {
  const out = [];
  for (const ent of readdirSync(dir, { withFileTypes: true })) {
    const full = join(dir, ent.name);
    if (ent.isDirectory()) {
      out.push(...listSourceFiles(full));
    } else if (/\.(ts|tsx)$/.test(ent.name)) {
      out.push(full);
    }
  }
  return out;
}

// POSIX-normalized path relative to src/ (so allowlist + exclusion are OS-independent).
function relToSrc(absPath) {
  return relative(SRC_DIR, absPath).split(sep).join("/");
}

// src/components/ui/** is excluded entirely (it is where appearance is authored).
function isExcluded(relPath) {
  return relPath.startsWith("components/ui/");
}

// Scan one file's text. Returns [{ rule, file, line, text }].
// Exported so the fixture test can drive it without touching disk.
export function scanContent(relPath, content) {
  if (isExcluded(relPath)) return [];
  const violations = [];
  const lines = content.split("\n");
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    for (const rule of RULES) {
      if (rule.allowlisted && rule.allowlisted(relPath)) continue;
      const matched = rule.test(line);
      if (matched != null) {
        violations.push({ rule: rule.id, file: relPath, line: i + 1, text: matched });
      }
    }
  }
  return violations;
}

// Content-hash key: (rule, normalized-violation-text) — NOT (file,line), NOT (file,rule)-count.
// Normalize: collapse all whitespace runs to one space + trim, so reflow/indent moves are no-ops.
export function hashViolation(rule, text) {
  const norm = text.replace(/\s+/g, " ").trim();
  return createHash("sha1").update(rule + " " + norm).digest("hex").slice(0, 16);
}

export function loadBaseline(path = BASELINE_PATH) {
  if (!existsSync(path)) return { hashes: [] };
  return JSON.parse(readFileSync(path, "utf8"));
}

// baseline_new ⊆ baseline_old check + not-in-baseline detection.
// Returns { notInBaseline: [...violations], grew: bool }.
export function diffBaseline(violations, baseline) {
  const baseSet = new Set(baseline.hashes ?? []);
  const seen = new Set();
  const notInBaseline = [];
  for (const v of violations) {
    const h = hashViolation(v.rule, v.text);
    seen.add(h);
    if (!baseSet.has(h)) notInBaseline.push({ ...v, hash: h });
  }
  // grew = the current set introduced a hash the baseline does not contain.
  return { notInBaseline, grew: notInBaseline.length > 0, seenHashes: seen };
}

// --- CLI ---------------------------------------------------------------------
function runCli(argv) {
  const initMode = argv.includes("--init");
  const files = listSourceFiles(SRC_DIR);
  const all = [];
  for (const abs of files) {
    const rel = relToSrc(abs);
    all.push(...scanContent(rel, readFileSync(abs, "utf8")));
  }

  if (initMode) {
    const hashes = [...new Set(all.map((v) => hashViolation(v.rule, v.text)))].sort();
    const baseline = {
      note: "Content-hash baseline for scripts/check-design.mjs. Keyed on (rule, normalized text). Ratchet: only ever shrinks. Regenerate seed with `node scripts/check-design.mjs --init`.",
      generated: new Date().toISOString().slice(0, 10),
      hashes,
    };
    writeFileSync(BASELINE_PATH, JSON.stringify(baseline, null, 2) + "\n");
    process.stderr.write(`[check-design] wrote ${hashes.length} baseline hashes to ${BASELINE_PATH}\n`);
    return 0;
  }

  const baseline = loadBaseline();
  const { notInBaseline } = diffBaseline(all, baseline);
  if (notInBaseline.length > 0) {
    process.stderr.write(`[check-design] ${notInBaseline.length} NEW design violation(s) not in baseline:\n`);
    for (const v of notInBaseline) {
      process.stderr.write(`  ${v.rule}  src/${v.file}:${v.line}  ${JSON.stringify(v.text)}  [${v.hash}]\n`);
    }
    process.stderr.write("Fix the appearance utility (move it into src/components/ui/**) or, if intentional and pre-existing, the baseline is wrong.\n");
    return 2;
  }
  return 0;
}

if (import.meta.url === `file://${process.argv[1]}`) {
  process.exit(runCli(process.argv.slice(2)));
}
