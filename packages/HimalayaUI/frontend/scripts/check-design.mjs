#!/usr/bin/env node
// Zero-dependency design-token guard. Globs src/**/*.{ts,tsx}, excludes src/components/ui/**,
// flags banned appearance utilities / raw color literals. PURE-ABSOLUTE: every matched
// violation (after the ui/ exclusion + per-rule allowlist) is a hard error (exit 2). There is
// no baseline — the named scale + component library are now the only sanctioned source of
// appearance, so any new bracket/raw-color escape is a regression, full stop.
// Pure functions are exported for unit testing; the CLI runs only when invoked directly.
import { readdirSync, readFileSync } from "node:fs";
import { join, relative, sep, posix } from "node:path";
import { fileURLToPath } from "node:url";

// import.meta.url is a file: URL when run by node directly, but some loaders
// (e.g. Vitest's transform) hand back a non-file scheme that fileURLToPath rejects.
// Fall back to a cwd-relative scripts dir so module load never throws; the CLI
// (which actually walks the tree) is only ever invoked via `node`, where the
// file: URL path holds.
function scriptsDir() {
  try {
    return fileURLToPath(new URL(".", import.meta.url));
  } catch {
    return join(process.cwd(), "scripts");
  }
}
const HERE = scriptsDir();
const SRC_DIR = join(HERE, "..", "src");

// --- Ban rules ---------------------------------------------------------------
// Each rule: { id, test(line) -> matched-substring|null }. Rules #3 and #5 also consult
// `relPath` (POSIX-normalized, relative to src/) against the shared color-authoring allowlist.
//
// Color-authoring allowlist: the files where color is legitimately authored (palette tables,
// computed coloring, figure-export presets, the canvas/detector layers that paint pixels, the
// app entry that seeds CSS vars). Rule #3 (arbitrary-value color UTILITIES that need a runtime
// color) already exempts them; rule #5 (raw color literals anywhere) MUST exempt the same set —
// otherwise the palette layer can never be clean, since authoring a palette necessarily writes
// oklch()/rgba() literals. (Rules #1/#2/#4 take NO allowlist.)
const COLOR_AUTHORING_ALLOWLIST = new Set([
  "phases.ts",
  "lib/comparison/coloring.ts",
  "components/MemberHeatmapLayer.tsx",
  "components/DetectorImage.tsx",
  "components/FocusDetectorPanel.tsx",
  "main.tsx",
]);
// figure-export/** is allowlisted by prefix (the whole export palette dir).
const COLOR_AUTHORING_ALLOW_PREFIXES = ["lib/figure-export/"];

function colorAuthoringAllowed(relPath) {
  if (COLOR_AUTHORING_ALLOWLIST.has(relPath)) return true;
  return COLOR_AUTHORING_ALLOW_PREFIXES.some((p) => relPath.startsWith(p));
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
    // shadow-[…] is excluded — legit Plate-Lift rgba shadows live in tsx.
    id: "no-raw-color-utility",
    allowlisted: colorAuthoringAllowed,
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
    // Strip every var(--color-*) token AND every shadow-[…] value FIRST, then flag a
    // remaining raw literal. Stripping var() passes legit computed-color sites; stripping
    // shadow-[…] keeps rule #5 consistent with rule #3, which deliberately exempts the
    // Plate-Lift rgba shadows (the rgba lives inside an elevation token, not a color role).
    // What remains and trips: the raw oklch(0.05…) scrim hand-inlines and quoted hex colors.
    id: "no-raw-color-literal",
    allowlisted: colorAuthoringAllowed,
    test: (line) => {
      const stripped = line
        .replace(/var\(--color-[^)]*\)/g, "")
        .replace(/shadow-\[[^\]]*\]/g, "");
      const re = /(?:oklch|oklab|hsla?|rgba?)\(|["'`]#[0-9a-fA-F]{3,8}\b/;
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

// src/components/ui/** and src/print/ui/** are excluded entirely (where appearance is authored).
function isExcluded(relPath) {
  return relPath.startsWith("components/ui/") || relPath.startsWith("print/ui/");
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

// Relative import/export specifiers on one line. Group 1 = `... from "X"`, group 2 = dynamic import("X").
const IMPORT_SPEC_RE =
  /(?:import|export)\b[^'"]*?\bfrom\s*["']([^"']+)["']|import\s*\(\s*["']([^"']+)["']\s*\)/g;

// Import-boundary guard: a file under src/print/** may not import (relatively) from the OLD
// src/components/** or src/pages/**. We resolve each relative specifier against the importer's
// dir (POSIX) and flag it if it lands under components/ or pages/. Print-internal imports
// (./components, ./pages — i.e. src/print/components, src/print/pages) resolve under print/ and pass.
// Exported pure for unit testing.
export function scanLegacyImports(relPath, content) {
  if (!relPath.startsWith("print/")) return [];
  const dir = posix.dirname(relPath);
  const out = [];
  const lines = content.split("\n");
  for (let i = 0; i < lines.length; i++) {
    IMPORT_SPEC_RE.lastIndex = 0;
    let m;
    while ((m = IMPORT_SPEC_RE.exec(lines[i])) != null) {
      const spec = m[1] ?? m[2];
      if (!spec || !spec.startsWith(".")) continue; // bare/external specifiers are fine
      const resolved = posix.normalize(posix.join(dir, spec));
      if (resolved.startsWith("components/") || resolved.startsWith("pages/")) {
        out.push({ rule: "no-legacy-import", file: relPath, line: i + 1, text: spec });
      }
    }
  }
  return out;
}

// --- CLI ---------------------------------------------------------------------
// Pure-absolute: scan every source file, error on ANY violation (no baseline).
function runCli() {
  const files = listSourceFiles(SRC_DIR);
  const all = [];
  for (const abs of files) {
    const rel = relToSrc(abs);
    const content = readFileSync(abs, "utf8");
    all.push(...scanContent(rel, content));
    all.push(...scanLegacyImports(rel, content));
  }

  if (all.length > 0) {
    process.stderr.write(`[check-design] ${all.length} design violation(s):\n`);
    for (const v of all) {
      process.stderr.write(`  ${v.rule}  src/${v.file}:${v.line}  ${JSON.stringify(v.text)}\n`);
    }
    process.stderr.write(
      "Move the appearance utility into src/components/ui/** (or src/print/ui/**), or use a named " +
        "scale/role token. Raw color literals belong only in the color-authoring files. " +
        "src/print/** may not import from src/components/** or src/pages/** (no-legacy-import).\n",
    );
    return 2;
  }
  return 0;
}

if (import.meta.url === `file://${process.argv[1]}`) {
  process.exit(runCli());
}
