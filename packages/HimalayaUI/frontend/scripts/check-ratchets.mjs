#!/usr/bin/env node
import { readFileSync, writeFileSync, readdirSync } from "node:fs";
import { join, relative, sep, dirname } from "node:path";
import { fileURLToPath } from "node:url";

const SCRIPTS_DIR = (() => {
  try {
    return dirname(fileURLToPath(import.meta.url));
  } catch {
    return join(process.cwd(), "scripts");
  }
})();
const PKG_DIR = join(SCRIPTS_DIR, "..", "..");
const FRONTEND_SRC = join(PKG_DIR, "frontend", "src");
const BACKEND_SRC = join(PKG_DIR, "src");
const BASELINE_PATH = join(SCRIPTS_DIR, "ratchet-baseline.json");

const APPEARANCE_EXEMPT_PREFIXES = ["print/ui/", "print/plot/", "print/detector/", "print/comb/", "print/export/"];

const APPEARANCE_TOKEN =
  /(?<![\w-])(text-(?!left\b|right\b|center\b|justify\b|start\b|end\b|wrap\b|nowrap\b|balance\b|pretty\b|ellipsis\b|clip\b)[a-z][\w-]*|bg-[\w-]+|border(?:-[\w-]+)?|rounded(?:-[\w-]+)?|shadow(?:-[\w-]+)?|font-[\w-]+|tracking-[\w-]+|leading-[\w-]+|uppercase|lowercase|italic|underline|opacity-\d+|ring(?:-[\w-]+)?|(?:hover|focus|focus-visible|active):[\w-]+)(?![\w-])/;

const CLASSNAME_STRING = /className=(?:"([^"]*)"|\{`([^`]*)`\}|\{[^}]*?"([^"]*)")/g;

const COMMENT_EXEMPT = /eslint-disable|@ts-expect-error|@ts-ignore|physics:/;

export function countAppearance(relPath, content) {
  if (relPath.includes(".stories.")) return 0;
  if (APPEARANCE_EXEMPT_PREFIXES.some((p) => relPath.startsWith(p))) return 0;
  let n = 0;
  for (const m of content.matchAll(CLASSNAME_STRING)) {
    const cls = m[1] ?? m[2] ?? m[3] ?? "";
    if (APPEARANCE_TOKEN.test(cls)) n++;
  }
  return n;
}

export function countTsComments(content) {
  let n = 0;
  let inBlock = false;
  for (const raw of content.split("\n")) {
    const line = raw.trim();
    if (inBlock) {
      if (!COMMENT_EXEMPT.test(line)) n++;
      if (line.includes("*/")) inBlock = false;
      continue;
    }
    const isLine = line.startsWith("//");
    const isBlock = line.startsWith("/*") || line.startsWith("{/*");
    if (!isLine && !isBlock) continue;
    if (!COMMENT_EXEMPT.test(line)) n++;
    if (isBlock && !line.includes("*/")) inBlock = true;
  }
  return n;
}

export function countJuliaComments(content) {
  let n = 0;
  let inDoc = false;
  for (const raw of content.split("\n")) {
    if (inDoc) {
      n++;
      if (raw.trim().endsWith('"""')) inDoc = false;
      continue;
    }
    if (raw.startsWith('"""')) {
      n++;
      const rest = raw.slice(3);
      if (!rest.includes('"""')) inDoc = true;
      continue;
    }
    const line = raw.trim();
    if (line.startsWith("#") && !COMMENT_EXEMPT.test(line)) n++;
  }
  return n;
}

export function compareToBaseline(baseline, current) {
  const increased = [];
  const decreased = [];
  for (const kind of ["appearance", "comments"]) {
    const base = baseline[kind] ?? {};
    const cur = current[kind] ?? {};
    const files = [...new Set([...Object.keys(base), ...Object.keys(cur)])].sort();
    for (const file of files) {
      const was = base[file] ?? 0;
      const now = cur[file] ?? 0;
      if (now > was) increased.push({ kind, file, was, now });
      else if (now < was) decreased.push({ kind, file, was, now });
    }
  }
  return { increased, decreased };
}

function walk(dir, exts, out = []) {
  for (const ent of readdirSync(dir, { withFileTypes: true })) {
    const p = join(dir, ent.name);
    if (ent.isDirectory()) walk(p, exts, out);
    else if (exts.some((e) => ent.name.endsWith(e))) out.push(p);
  }
  return out;
}

const posix = (p) => p.split(sep).join("/");

export function measure() {
  const appearance = {};
  const comments = {};
  for (const abs of walk(FRONTEND_SRC, [".ts", ".tsx"])) {
    const rel = posix(relative(FRONTEND_SRC, abs));
    if (rel.includes(".stories.")) continue;
    const key = `frontend/src/${rel}`;
    const content = readFileSync(abs, "utf8");
    const a = countAppearance(rel, content);
    const c = countTsComments(content);
    if (a > 0) appearance[key] = a;
    if (c > 0) comments[key] = c;
  }
  for (const abs of walk(BACKEND_SRC, [".jl"])) {
    const key = `src/${posix(relative(BACKEND_SRC, abs))}`;
    const c = countJuliaComments(readFileSync(abs, "utf8"));
    if (c > 0) comments[key] = c;
  }
  return { appearance, comments };
}

const total = (m) => Object.values(m).reduce((s, v) => s + v, 0);

const FIX = {
  appearance:
    "Styling belongs in a print/ui primitive. Use an existing primitive or add a variant to one; a consumer className carries placement only (margin, position, grid).",
  comments:
    "HimalayaUI code is comment-free. Put the reason in a name, a test or an assertion. Only `// physics:` / `# physics:` lines and tool directives (eslint-disable with a reason, @ts-expect-error) are allowed.",
};

function main() {
  const update = process.argv.includes("--update");
  const baseline = JSON.parse(readFileSync(BASELINE_PATH, "utf8"));
  const current = measure();
  const bootstrap = Object.keys(baseline).length === 0;
  const { increased, decreased } = compareToBaseline(baseline, current);

  if (update && bootstrap) {
    writeFileSync(BASELINE_PATH, `${JSON.stringify(current, null, 2)}\n`);
    console.warn(
      `ratchet: baseline created (appearance ${total(current.appearance)}, comments ${total(current.comments)})`,
    );
    return;
  }

  if (increased.length > 0) {
    for (const kind of ["appearance", "comments"]) {
      const rows = increased.filter((r) => r.kind === kind);
      if (rows.length === 0) continue;
      console.error(`\nratchet: ${kind} grew in ${rows.length} file(s):`);
      for (const r of rows) console.error(`  ${r.file}: ${r.was} -> ${r.now}`);
      console.error(`  Fix: ${FIX[kind]}`);
    }
    process.exit(2);
  }

  if (update) {
    writeFileSync(BASELINE_PATH, `${JSON.stringify(current, null, 2)}\n`);
    console.warn(
      `ratchet: baseline updated (appearance ${total(current.appearance)}, comments ${total(current.comments)})`,
    );
    return;
  }

  if (decreased.length > 0) {
    console.error(`\nratchet: counts fell in ${decreased.length} file(s). Lock in the progress:`);
    console.error("  npm run lint:ratchet -- --update   (then commit scripts/ratchet-baseline.json)");
    process.exit(2);
  }

  console.warn(
    `ratchet: ok (appearance ${total(current.appearance)}, comments ${total(current.comments)})`,
  );
}

const invokedDirectly = (() => {
  try {
    return process.argv[1] !== undefined && fileURLToPath(import.meta.url) === process.argv[1];
  } catch {
    return false;
  }
})();
if (invokedDirectly) main();
