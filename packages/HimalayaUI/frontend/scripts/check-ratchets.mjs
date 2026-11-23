#!/usr/bin/env node
import { readFileSync, writeFileSync, existsSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { join, relative, sep, dirname } from "node:path";
import { fileURLToPath } from "node:url";
import ts from "typescript";

const SCRIPTS_DIR = dirname(fileURLToPath(import.meta.url));
const PKG_DIR = join(SCRIPTS_DIR, "..", "..");
const FRONTEND_SRC = join(PKG_DIR, "frontend", "src");
const REPO = join(PKG_DIR, "..", "..");

export const TS_EXTENSIONS = [".ts", ".tsx", ".mts", ".cts", ".js", ".jsx", ".mjs", ".cjs"];
export const JULIA_EXTENSIONS = [".jl"];
export const CSS_EXTENSIONS = [".css"];
export const EXCLUDED = [];
const BASELINE_PATH = join(SCRIPTS_DIR, "ratchet-baseline.json");

const APPEARANCE_EXEMPT_PREFIXES = ["ui/", "plot/", "detector/"];

const APPEARANCE_TOKEN =
  /^(text-(?!(?:left|right|center|justify|start|end|wrap|nowrap|balance|pretty|ellipsis|clip)$)[a-z][\w-]*|bg-[\w-]+|border(?:-[\w-]+)?|rounded(?:-[\w-]+)?|shadow(?:-[\w-]+)?|font-[\w-]+|tracking-[\w-]+|leading-[\w-]+|uppercase|lowercase|italic|underline|opacity-\d+|ring(?:-[\w-]+)?)$/;

const TS_EXEMPT = /^(?:\/\/\/\s*<reference\b|(?:\/\/+|\/\*+)\s*(?:eslint-disable[\w-]*\s+\S.*?\s--\s+\S|@ts-expect-error\b|physics:))/s;
const DIRECTIVE = /^(?:\/\/\/|\/\/+|\/\*+)\s*(?:<reference\b|eslint-disable|@ts-expect-error)/;
export const DIRECTIVE_MAX = 120;

function isExemptTsComment(text) {
  if (!TS_EXEMPT.test(text)) return false;
  if (!DIRECTIVE.test(text)) return true;
  return !text.includes("\n") && text.length <= DIRECTIVE_MAX;
}
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

const isAppearanceToken = (tok) =>
  APPEARANCE_TOKEN.test((tok.split(":").pop() ?? "").replace(/^!/, "").replace(/\/[\w.[\]]+$/, ""));

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
    if (isExemptTsComment(content.slice(pos, end))) return;
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
const DEFINITION =
  /^[ \t]*(?:function\b|macro\b|struct\b|mutable\s+struct\b|abstract\s+type\b|primitive\s+type\b|const\b|module\b|baremodule\b|global\b|@\w|[^\s=#][^\n=]*\([^\n]*\)[^\n=]*=(?!=)|[A-Za-z_][\w!]*\s*::|(?!(?:end|else|elseif|catch|finally|return|nothing|true|false|break|continue)\b)[A-Za-z_][\w!]*[ \t]*(?:\n|$))/;
const IDENT_START = /[A-Za-z_]/;
const IDENT = /[\w!]/;

function skipJlComment(s, i) {
  if (s[i + 1] !== "=") {
    while (i < s.length && s[i] !== "\n") i++;
    return i;
  }
  let depth = 1;
  i += 2;
  while (i < s.length && depth > 0) {
    if (s.startsWith("#=", i)) {
      depth++;
      i += 2;
    } else if (s.startsWith("=#", i)) {
      depth--;
      i += 2;
    } else i++;
  }
  return i;
}

function skipJlCode(s, i, close) {
  let depth = 0;
  while (i < s.length) {
    const c = s[i];
    if (c === '"' || c === "`") i = skipJlDelimited(s, i);
    else if (c === "#") i = skipJlComment(s, i);
    else if (c === "(" || c === "[" || c === "{") {
      depth++;
      i++;
    } else if (c === ")" || c === "]" || c === "}") {
      if (depth === 0 && c === close) return i + 1;
      depth--;
      i++;
    } else i++;
  }
  return i;
}

function skipJlDelimited(s, i) {
  const q = s[i];
  const fence = s.startsWith(q.repeat(3), i) ? q.repeat(3) : q;
  i += fence.length;
  while (i < s.length && !s.startsWith(fence, i)) {
    if (s[i] === "\\") i += 2;
    else if (s[i] === "$" && s[i + 1] === "(") i = skipJlCode(s, i + 2, ")");
    else i++;
  }
  return i + fence.length;
}

function restOfLineIsBlank(s, i) {
  while (i < s.length && (s[i] === " " || s[i] === "\t" || s[i] === "\r")) i++;
  return i >= s.length || s[i] === "\n" || s[i] === "#";
}

function nextLineDefines(s, i) {
  while (i < s.length && s[i] !== "\n") i++;
  while (i < s.length && /\s/.test(s[i])) i++;
  let j = i;
  while (j > 0 && s[j - 1] !== "\n") j--;
  return i < s.length && DEFINITION.test(s.slice(j, s.indexOf("\n", i) === -1 ? s.length : s.indexOf("\n", i) + 1));
}

export function countJuliaComments(content) {
  const s = content;
  const starts = lineStarts(s);
  const lines = new Set();
  const mark = (from, to) => {
    for (let l = lineOf(starts, from); l <= lineOf(starts, Math.max(from, to - 1)); l++) lines.add(l);
  };
  let i = 0;
  let last = "";
  let atLineStart = true;
  while (i < s.length) {
    const c = s[i];
    if (c === "\n") {
      atLineStart = true;
      i++;
      continue;
    }
    if (c === " " || c === "\t" || c === "\r") {
      i++;
      continue;
    }
    const statementStart = atLineStart && !CONTINUES.has(last);
    atLineStart = false;
    if (c === "#") {
      const end = skipJlComment(s, i);
      if (!JL_EXEMPT.test(s.slice(i, end))) mark(i, end);
      i = end;
      continue;
    }
    const docMacro = /^@doc\s+(?=(?:raw)?")/.exec(s.slice(i, i + 16));
    const rawPrefix = /^raw(?=")/.exec(s.slice(i, i + 4));
    if (docMacro || rawPrefix || c === '"') {
      const quote = i + (docMacro ? docMacro[0].length : 0) + (rawPrefix ? 3 : 0);
      const q = s.startsWith("raw", quote) ? quote + 3 : quote;
      const end = skipJlDelimited(s, q);
      const triple = s.startsWith('"""', q);
      const isDoc =
        statementStart &&
        (docMacro !== null || triple || (restOfLineIsBlank(s, end) && nextLineDefines(s, end)));
      if (isDoc) mark(i, end);
      i = end;
      last = '"';
      continue;
    }
    if (c === "`") {
      i = skipJlDelimited(s, i);
      last = "`";
      continue;
    }
    if (c === "'" && !/[\w)\]'.]/.test(last)) {
      i++;
      while (i < s.length && s[i] !== "'" && s[i] !== "\n") i += s[i] === "\\" ? 2 : 1;
      i++;
      last = "'";
      continue;
    }
    if (IDENT_START.test(c)) {
      let j = i + 1;
      while (j < s.length && IDENT.test(s[j])) j++;
      if (s[j] === '"' || s[j] === "`") {
        i = skipJlDelimited(s, j);
        last = '"';
      } else {
        last = s[j - 1];
        i = j;
      }
      continue;
    }
    last = c;
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

function renamesCover(base, cur, exists) {
  const vanished = Object.entries(base)
    .filter(([file]) => !(file in cur) && !exists(file))
    .map(([, n]) => n);
  for (const [file, n] of Object.entries(cur)) {
    if ((base[file] ?? 0) >= n) continue;
    if (file in base) return false;
    const k = vanished.indexOf(n);
    if (k === -1) return false;
    vanished.splice(k, 1);
  }
  return true;
}

export function canUpdate(baseline, current, exists = (file) => existsSync(join(REPO, file))) {
  return ["appearance", "comments"].every((k) => renamesCover(baseline[k] ?? {}, current[k] ?? {}, exists));
}

const posix = (p) => p.split(sep).join("/");

export function countCssComments(content) {
  const starts = lineStarts(content);
  const lines = new Set();
  let i = 0;
  while (i < content.length) {
    const c = content[i];
    if (c === '"' || c === "'") {
      i++;
      while (i < content.length && content[i] !== c && content[i] !== "\n") i += content[i] === "\\" ? 2 : 1;
      i++;
    } else if (content.startsWith("/*", i)) {
      const close = content.indexOf("*/", i + 2);
      const end = close === -1 ? content.length : close + 2;
      if (!/^\/\*+\s*physics:/.test(content.slice(i, end))) {
        for (let l = lineOf(starts, i); l <= lineOf(starts, end - 1); l++) lines.add(l);
      }
      i = end;
    } else i++;
  }
  return lines.size;
}

const hasExtension = (file, exts) => exts.some((e) => file.endsWith(e));

export function trackedSources() {
  const out = execFileSync("git", ["ls-files", "-z"], { cwd: REPO, encoding: "utf8", maxBuffer: 64 * 1024 * 1024 });
  return out
    .split("\0")
    .filter((f) => f !== "" && hasExtension(f, [...TS_EXTENSIONS, ...JULIA_EXTENSIONS, ...CSS_EXTENSIONS]))
    .filter((f) => !EXCLUDED.includes(f))
    .filter((f) => existsSync(join(REPO, f)));
}

function countComments(file, content) {
  if (hasExtension(file, JULIA_EXTENSIONS)) return countJuliaComments(content);
  if (hasExtension(file, CSS_EXTENSIONS)) return countCssComments(content);
  return countTsComments(content);
}

export function measure() {
  const appearance = {};
  const comments = {};
  const frontendSrc = posix(relative(REPO, FRONTEND_SRC)) + "/";
  for (const file of trackedSources()) {
    const content = readFileSync(join(REPO, file), "utf8");
    const c = countComments(file, content);
    if (c > 0) comments[file] = c;
    if (file.startsWith(frontendSrc) && hasExtension(file, [".ts", ".tsx"])) {
      const rel = file.slice(frontendSrc.length);
      const a = countAppearance(rel, content);
      if (a > 0) appearance[`frontend/src/${rel}`] = a;
    }
  }
  return { appearance, comments };
}

const FIX = {
  appearance:
    "Styling belongs in a src/ui primitive. Use an existing primitive or add a variant to one; a consumer className carries placement only (margin, position, grid).",
  comments:
    "Code in this repo is comment-free, tests, tools and stories included. Put the reason in a clear name, a small function, a type, a test or an assertion. Allowed: comments that start with `physics:`, and one-line `eslint-disable… -- <reason>`, `@ts-expect-error` and `/// <reference>` directives of at most 120 characters.",
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
      console.error(`\nratchet: ${kind} grew in ${rows.length} file(s):`);
      for (const r of rows) console.error(`  ${r.file}: ${r.was} -> ${r.now}`);
      if (canUpdate(baseline, current)) console.error("  These are renames of files with the same count. Run `npm run lint:ratchet -- --update`.");
      else console.error(`  Fix: ${FIX[kind]} A split or rename that keeps comments counts as growth: delete them.`);
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
