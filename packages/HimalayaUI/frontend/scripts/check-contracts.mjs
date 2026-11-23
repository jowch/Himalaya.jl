#!/usr/bin/env node
import { readFileSync, existsSync } from "node:fs";
import { join, dirname } from "node:path";
import { fileURLToPath } from "node:url";

const FRONTEND = join(dirname(fileURLToPath(import.meta.url)), "..");
const REPO = join(FRONTEND, "..", "..", "..");
const STATUSES = new Set(["built", "untested", "superseded", "planned"]);
const ID = /^[a-z0-9]+(?:-[a-z0-9]+)*$/;
const DATE = /^\d{4}-\d{2}-\d{2}$/;

export const MOCKUP = { name: "docs/mockup-contract.md", pinned: true, dated: false, root: FRONTEND };
export const DECISIONS = { name: "docs/decisions.md", pinned: false, dated: true, root: REPO };

function cells(line) {
  return line.trim().replace(/^\|/, "").replace(/\|$/, "").split("|").map((c) => c.trim());
}

export function parseContract(text) {
  const lines = text.split("\n");
  const version = /Mock-up version: `([^`]+)`/.exec(text)?.[1];
  const sha256 = /Mock-up sha256: `([0-9a-f]{64})`/.exec(text)?.[1];
  const rows = [];
  let columns = null;
  lines.forEach((line, i) => {
    if (!line.trim().startsWith("|")) {
      columns = null;
      return;
    }
    const row = cells(line);
    if (row[0] === "id") {
      columns = Object.fromEntries(row.map((c, k) => [c.toLowerCase(), k]));
      return;
    }
    const [id] = row;
    if (!id || !ID.test(id)) return;
    const at = (name, fallback) => row[columns?.[name] ?? fallback] ?? "";
    const parsed = { id, behaviour: row[1] ?? "", status: at("status", 2), evidence: at("evidence", 3), line: i + 1 };
    if (columns && "decided" in columns) parsed.decided = at("decided");
    rows.push(parsed);
  });
  return { version, sha256, rows };
}

const SKIPPED = new Set(["skip", "todo", "only", "fixme", "fails"]);

function escape(s) {
  return s.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

function declaration(content, title) {
  const re = new RegExp(`\\b(?:it|test)(?:\\.(\\w+))?\\s*\\(\\s*(["'\`])${escape(title)}\\2`, "g");
  let found;
  for (const m of content.matchAll(re)) {
    if (!found || (found.modifier && !m[1])) found = { modifier: m[1] };
  }
  if (!found && new RegExp(`@testset\\s+"${escape(title)}"`).test(content)) found = {};
  return found;
}

function builtProblem(row, read) {
  const m = /^(\S+) › (.+)$/.exec(row.evidence);
  if (!m) return "a built row's evidence must be `path › test title`";
  const [, path, title] = m;
  const content = read(path);
  if (content === undefined) return `${path} does not exist`;
  const decl = declaration(content, title);
  if (!decl) return `${path} does not contain the test "${title}"`;
  if (decl.modifier && SKIPPED.has(decl.modifier)) {
    return `the test "${title}" is marked .${decl.modifier}, so it does not run as written; cite a test that runs`;
  }
  return null;
}

export function checkContract(contract, read, doc = MOCKUP) {
  const problems = [];
  if (doc.pinned && !contract.version) problems.push(`${doc.name} does not pin the mock-up version`);
  if (doc.pinned && !contract.sha256) problems.push(`${doc.name} does not pin the mock-up sha256`);
  const seen = new Set();
  for (const row of contract.rows) {
    const at = `${doc.name}:${row.line} ${row.id}`;
    if (seen.has(row.id)) problems.push(`${at}: the id is used more than once`);
    seen.add(row.id);
    if (doc.dated && !DATE.test(row.decided ?? "")) problems.push(`${at}: the Decided column must be a YYYY-MM-DD date`);
    if (!STATUSES.has(row.status)) {
      problems.push(`${at}: status "${row.status}" is not one of ${[...STATUSES].join(", ")}`);
      continue;
    }
    if (row.evidence === "") {
      problems.push(`${at}: no evidence`);
      continue;
    }
    if (row.status === "built") {
      const p = builtProblem(row, read);
      if (p) problems.push(`${at}: ${p}`);
    } else if (row.status === "superseded" && !/\b\d{4}-\d{2}-\d{2}\b|\bPR \d+|docs\//.test(row.evidence)) {
      problems.push(`${at}: a superseded row must cite its decision with a date, a PR number or a doc path`);
    } else if (row.status === "planned" && !/\bPR \d+|#\d+|\bthread "/.test(row.evidence)) {
      problems.push(`${at}: a planned row must name its owner: a PR number, an issue or a thread`);
    }
  }
  return problems;
}

function reader(root) {
  return (path) => {
    const full = join(root, path);
    return existsSync(full) ? readFileSync(full, "utf8") : undefined;
  };
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  let failed = 0;
  for (const doc of [MOCKUP, DECISIONS]) {
    const contract = parseContract(readFileSync(join(REPO, doc.name), "utf8"));
    const problems = checkContract(contract, reader(doc.root), doc);
    const counts = {};
    for (const r of contract.rows) counts[r.status] = (counts[r.status] ?? 0) + 1;
    if (problems.length) {
      for (const p of problems) console.error(p);
      console.error(`\n${problems.length} problem(s) in ${doc.name}. Each row needs a status and the evidence that status requires; see the Statuses table at the top of the file.`);
      failed += problems.length;
      continue;
    }
    const pin = doc.pinned ? ` against mock-up ${contract.version}` : "";
    console.log(`${doc.name}: ${contract.rows.length} rows (${Object.entries(counts).map(([k, v]) => `${v} ${k}`).join(", ")})${pin}`);
  }
  if (failed) process.exit(1);
}
