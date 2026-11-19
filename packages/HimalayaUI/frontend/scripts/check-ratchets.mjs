#!/usr/bin/env node
import { readFileSync, writeFileSync, readdirSync } from "node:fs";
import { join, relative, sep, dirname } from "node:path";
import { fileURLToPath } from "node:url";
import ts from "typescript";

const SCRIPTS_DIR = dirname(fileURLToPath(import.meta.url));
const PKG_DIR = join(SCRIPTS_DIR, "..", "..");
const FRONTEND_SRC = join(PKG_DIR, "frontend", "src");
const BACKEND_SRC = join(PKG_DIR, "src");
const BASELINE_PATH = join(SCRIPTS_DIR, "ratchet-baseline.json");

const APPEARANCE_EXEMPT_PREFIXES = ["print/ui/", "print/plot/", "print/detector/", "print/comb/", "print/export/"];

const APPEARANCE_TOKEN =
  /^(text-(?!(?:left|right|center|justify|start|end|wrap|nowrap|balance|pretty|ellipsis|clip)$)[a-z][\w-]*|bg-[\w-]+|border(?:-[\w-]+)?|rounded(?:-[\w-]+)?|shadow(?:-[\w-]+)?|font-[\w-]+|tracking-[\w-]+|leading-[\w-]+|uppercase|lowercase|italic|underline|opacity-\d+|ring(?:-[\w-]+)?)$/;

const TS_EXEMPT = /^(?:\/\/+|\/\*+)\s*(?:eslint-disable[\w-]*\s+\S.*?\s--\s+\S|@ts-expect-error\b|physics:)/s;
const JL_EXEMPT = /^#\s*physics:/;

function readQuoted(s, i) {
  const q = s[i];
  let j = i + 1;
  while (j < s.length && s[j] !== q) j += s[j] === "\\" ? 2 : 1;
  return { text: s.slice(i + 1, j), end: j + 1 };
}

function readTemplate(s, i, out) {
  let j = i + 1;
  let chunk = "";
  while (j < s.length && s[j] !== "`") {
    if (s[j] === "\\") {
      chunk += s.slice(j, j + 2);
      j += 2;
    } else if (s[j] === "$" && s[j + 1] === "{") {
      out.push(chunk);
      chunk = "";
      j = readExpression(s, j + 2, "}", out);
    } else {
      chunk += s[j];
      j++;
    }
  }
  out.push(chunk);
  return j + 1;
}

function readExpression(s, i, close, out) {
  let depth = 0;
  let j = i;
  while (j < s.length) {
    const c = s[j];
    if (c === '"' || c === "'") {
      const { text, end } = readQuoted(s, j);
      out.push(text);
      j = end;
    } else if (c === "`") {
      j = readTemplate(s, j, out);
    } else if (c === "{" || c === "(" || c === "[") {
      depth++;
      j++;
    } else if (c === "}" || c === ")" || c === "]") {
      if (depth === 0 && c === close) return j + 1;
      depth--;
      j++;
    } else {
      j++;
    }
  }
  return j;
}

export function classNameLiterals(content) {
  const out = [];
  const re = /className=/g;
  let m;
  while ((m = re.exec(content)) !== null) {
    const i = m.index + m[0].length;
    if (content[i] === '"' || content[i] === "'") {
      out.push(readQuoted(content, i).text);
    } else if (content[i] === "{") {
      readExpression(content, i + 1, "}", out);
    }
  }
  return out;
}

const isAppearanceToken = (tok) => APPEARANCE_TOKEN.test(tok.split(":").pop() ?? "");

export function countAppearance(relPath, content) {
  if (relPath.includes(".stories.")) return 0;
  if (APPEARANCE_EXEMPT_PREFIXES.some((p) => relPath.startsWith(p))) return 0;
  let n = 0;
  for (const lit of classNameLiterals(content)) {
    for (const tok of lit.split(/\s+/)) if (tok && isAppearanceToken(tok)) n++;
  }
  return n;
}

function lineStarts(content) {
  const starts = [0];
  for (let i = 0; i < content.length; i++) if (content[i] === "\n") starts.push(i + 1);
  return starts;
}

function lineOf(starts, pos) {
  let lo = 0;
  let hi = starts.length - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (starts[mid] <= pos) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

export function countTsComments(content) {
  const sf = ts.createSourceFile("x.tsx", content, ts.ScriptTarget.Latest, true, ts.ScriptKind.TSX);
  const seen = new Set();
  const lines = new Set();
  const starts = lineStarts(content);
  const record = (pos, end) => {
    if (seen.has(pos)) return;
    seen.add(pos);
    if (TS_EXEMPT.test(content.slice(pos, end))) return;
    for (let l = lineOf(starts, pos); l <= lineOf(starts, end - 1); l++) lines.add(l);
  };
  const visit = (node) => {
    ts.forEachLeadingCommentRange(content, node.pos, record);
    ts.forEachTrailingCommentRange(content, node.end, record);
    for (const child of node.getChildren(sf)) visit(child);
  };
  visit(sf);
  return lines.size;
}

const CONTINUES = new Set(["(", ",", "=", "[", "{", "+", "-", "*", "/", "|", "&", "<", ">", "~", ":", "?", "$"]);

export function countJuliaComments(content) {
  const lines = new Set();
  let line = 0;
  let i = 0;
  let last = "";
  let atLineStart = true;
  const n = content.length;
  const advance = () => {
    if (content[i] === "\n") {
      line++;
      atLineStart = true;
    }
    i++;
  };
  while (i < n) {
    const c = content[i];
    if (c === "\n") {
      advance();
      continue;
    }
    if (atLineStart && (c === " " || c === "\t")) {
      i++;
      continue;
    }
    const lineStart = atLineStart;
    atLineStart = false;
    if (c === "#" && content[i + 1] === "=") {
      let depth = 1;
      lines.add(line);
      i += 2;
      while (i < n && depth > 0) {
        if (content[i] === "#" && content[i + 1] === "=") {
          depth++;
          i += 2;
        } else if (content[i] === "=" && content[i + 1] === "#") {
          depth--;
          i += 2;
        } else {
          advance();
          atLineStart = false;
        }
        lines.add(line);
      }
      continue;
    }
    if (c === "#") {
      let j = i;
      while (j < n && content[j] !== "\n") j++;
      if (!JL_EXEMPT.test(content.slice(i, j))) lines.add(line);
      i = j;
      continue;
    }
    if (content.startsWith('"""', i)) {
      const isDoc = lineStart && !CONTINUES.has(last);
      const startLine = line;
      i += 3;
      while (i < n && !content.startsWith('"""', i)) {
        if (content[i] === "\\") i++;
        advance();
        atLineStart = false;
      }
      i += 3;
      if (isDoc) for (let l = startLine; l <= line; l++) lines.add(l);
      last = '"';
      continue;
    }
    if (c === '"') {
      i++;
      while (i < n && content[i] !== '"') {
        if (content[i] === "\\") i++;
        advance();
        atLineStart = false;
      }
      i++;
      last = '"';
      continue;
    }
    if (c === "'" && !/[\w)\]'.]/.test(last)) {
      i++;
      while (i < n && content[i] !== "'" && content[i] !== "\n") i += content[i] === "\\" ? 2 : 1;
      i++;
      last = "'";
      continue;
    }
    if (c !== " " && c !== "\t" && c !== "\r") last = c;
    i++;
  }
  return lines.size;
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

const total = (m) => Object.values(m ?? {}).reduce((s, v) => s + v, 0);

export function canUpdate(baseline, current) {
  return ["appearance", "comments"].every((k) => total(current[k]) <= total(baseline[k]));
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

const FIX = {
  appearance:
    "Styling belongs in a print/ui primitive. Use an existing primitive or add a variant to one; a consumer className carries placement only (margin, position, grid).",
  comments:
    "HimalayaUI code is comment-free. Put the reason in a name, a test or an assertion. Allowed: comments that start with `physics:`, `eslint-disable… -- <reason>`, and `@ts-expect-error`.",
};

function writeBaseline(current, verb) {
  writeFileSync(BASELINE_PATH, `${JSON.stringify(current, null, 2)}\n`);
  console.warn(
    `ratchet: baseline ${verb} (appearance ${total(current.appearance)} tokens, comments ${total(current.comments)} lines)`,
  );
}

function main() {
  const update = process.argv.includes("--update");
  const baseline = JSON.parse(readFileSync(BASELINE_PATH, "utf8"));
  const current = measure();

  if (update && Object.keys(baseline).length === 0) return writeBaseline(current, "created");
  if (update && canUpdate(baseline, current)) return writeBaseline(current, "updated");

  const { increased, decreased } = compareToBaseline(baseline, current);

  if (increased.length > 0) {
    for (const kind of ["appearance", "comments"]) {
      const rows = increased.filter((r) => r.kind === kind);
      if (rows.length === 0) continue;
      const grew = total(current[kind]) > total(baseline[kind]);
      console.error(`\nratchet: ${kind} grew in ${rows.length} file(s):`);
      for (const r of rows) console.error(`  ${r.file}: ${r.was} -> ${r.now}`);
      if (grew) console.error(`  Fix: ${FIX[kind]}`);
      else console.error("  The total did not grow (a rename or split). Run `npm run lint:ratchet -- --update`.");
    }
    process.exit(2);
  }

  if (decreased.length > 0) {
    console.error(`\nratchet: counts fell in ${decreased.length} file(s). Lock in the progress:`);
    console.error("  npm run lint:ratchet -- --update   (then commit scripts/ratchet-baseline.json)");
    process.exit(2);
  }

  console.warn(
    `ratchet: ok (appearance ${total(current.appearance)} tokens, comments ${total(current.comments)} lines)`,
  );
}

if (process.argv[1] !== undefined && fileURLToPath(import.meta.url) === process.argv[1]) main();
