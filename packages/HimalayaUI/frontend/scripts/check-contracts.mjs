#!/usr/bin/env node
import { readFileSync, existsSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { join, dirname } from "node:path";
import { fileURLToPath } from "node:url";

const FRONTEND = join(dirname(fileURLToPath(import.meta.url)), "..");
const REPO = join(FRONTEND, "..", "..", "..");
const STATUSES = new Set(["built", "untested", "superseded", "planned", "contradicted"]);
const ID = /^[a-z0-9]+(?:-[a-z0-9]+)*$/;
const DATE = /^\d{4}-\d{2}-\d{2}$/;

const LEADING_ID = /^(?:\d{4}-\d{2}-\d{2}[:,]?\s+)?([a-z][a-z0-9]*(?:-[a-z0-9]+)*)\b/;
const REF_ID = "([a-z][a-z0-9]*(?:-[a-z0-9]+)+)";
const DECISION_REFS = [new RegExp(`${REF_ID} in docs/decisions\\.md`, "g"), new RegExp(`docs/decisions\\.md ${REF_ID}`, "g")];

export const MOCKUP = { name: "docs/mockup-contract.md", pinned: true, dated: false, root: FRONTEND, floor: 70 };
export const DECISIONS = { name: "docs/decisions.md", pinned: false, dated: true, root: REPO, floor: 84, supersededById: true };

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

function blank(match) {
  return match.replace(/[^\n]/g, " ");
}

function withoutBlockComments(content) {
  return content.replace(/\/\*[\s\S]*?\*\//g, blank).replace(/#=[\s\S]*?=#/g, blank);
}

function declaration(content, title) {
  const code = withoutBlockComments(content);
  const re = new RegExp(`^[ \\t]*(?:it|test)(?:\\.(\\w+))?\\s*\\(\\s*(["'\`])${escape(title)}\\2`, "gm");
  let found;
  for (const m of code.matchAll(re)) {
    if (!found || (found.modifier && !m[1])) found = { modifier: m[1] };
  }
  if (!found && new RegExp(`^[ \\t]*@testset\\s+"${escape(title)}"`, "m").test(code)) found = {};
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

export function checkRowCount(contract, doc) {
  const n = contract.rows.length;
  if (n === 0) return [`${doc.name} has no rows. A merge or an edit has emptied it; restore it from main`];
  if (n < doc.floor) {
    return [`${doc.name} has ${n} rows, fewer than its floor of ${doc.floor}. A decision or behaviour that no longer applies is marked superseded, not deleted; if rows were removed on purpose, lower the floor in scripts/check-contracts.mjs in the same PR`];
  }
  return [];
}

export function decisionRefs(evidence) {
  return DECISION_REFS.flatMap((re) => [...evidence.matchAll(re)].map((m) => m[1]));
}

function supersededProblems(row, decisionIds) {
  if (!row.evidence.includes("docs/decisions.md")) return [];
  const refs = decisionRefs(row.evidence);
  if (refs.length === 0) return ["it mentions docs/decisions.md without naming a decision; cite it as `<id> in docs/decisions.md` or `docs/decisions.md <id>`"];
  return refs.filter((ref) => !decisionIds.has(ref)).map((ref) => `it cites ${ref} in docs/decisions.md, which has no row with that id`);
}

export function deletedRows(base, contract) {
  const now = new Set(contract.rows.map((r) => r.id));
  return base.rows.map((r) => r.id).filter((id) => !now.has(id));
}

export function checkContract(contract, read, doc = MOCKUP, decisionIds = undefined) {
  const problems = [];
  const ids = decisionIds ?? (doc.supersededById ? new Set(contract.rows.map((r) => r.id)) : undefined);
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
    } else if (row.status === "superseded" && doc.supersededById) {
      const by = LEADING_ID.exec(row.evidence)?.[1];
      if (!by) problems.push(`${at}: a superseded row starts its evidence with the id of the decision that replaced it, optionally after its date`);
      else if (by === row.id) problems.push(`${at}: a decision can't supersede itself`);
      else if (!ids.has(by)) problems.push(`${at}: it is superseded by ${by}, which has no row in ${doc.name}`);
    } else if (row.status === "superseded" && !/\b\d{4}-\d{2}-\d{2}\b|\bPR \d+|docs\//.test(row.evidence)) {
      problems.push(`${at}: a superseded row must cite its decision with a date, a PR number or a doc path`);
    } else if ((row.status === "planned" || row.status === "contradicted") && !/\bPR \d+|#\d+|\bthread "/.test(row.evidence)) {
      problems.push(`${at}: a ${row.status} row must name its owner: a PR number, an issue or a thread`);
    }
    if (row.status === "superseded" && ids) for (const p of supersededProblems(row, ids)) problems.push(`${at}: ${p}`);
  }
  return problems;
}

function reader(root) {
  return (path) => {
    const full = join(root, path);
    return existsSync(full) ? readFileSync(full, "utf8") : undefined;
  };
}

function git(args) {
  try {
    return execFileSync("git", args, { cwd: REPO, encoding: "utf8", stdio: ["ignore", "pipe", "ignore"] }).trim();
  } catch {
    return undefined;
  }
}

function baseRef() {
  const asked = process.env.CONTRACTS_BASE;
  if (asked) return git(["rev-parse", "--verify", `${asked}^{commit}`]) && asked;
  return git(["merge-base", "HEAD", "origin/main"]);
}

const BASE = process.argv[1] === fileURLToPath(import.meta.url) ? baseRef() : undefined;

function baseline(doc) {
  if (!BASE) return undefined;
  const text = git(["show", `${BASE}:${doc.name}`]);
  return text === undefined ? undefined : { ref: BASE.slice(0, 12), contract: parseContract(text) };
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  if (!BASE) console.log("No base commit to compare against (set CONTRACTS_BASE or fetch origin/main); checking row floors only.");
  let failed = 0;
  const load = (doc) => {
    const full = join(REPO, doc.name);
    return existsSync(full) ? parseContract(readFileSync(full, "utf8")) : undefined;
  };
  const decisions = load(DECISIONS);
  const decisionIds = decisions && new Set(decisions.rows.map((r) => r.id));
  for (const doc of [MOCKUP, DECISIONS]) {
    const contract = doc === DECISIONS ? decisions : load(doc);
    if (!contract) {
      console.error(`${doc.name} is missing. Restore it from main.`);
      failed += 1;
      continue;
    }
    const before = baseline(doc);
    const kept = before
      ? deletedRows(before.contract, contract).map((id) => `${doc.name}: the row ${id} is in ${before.ref} but not here. A row that no longer applies is marked superseded, not deleted`)
      : [];
    const problems = [...checkRowCount(contract, doc), ...kept, ...checkContract(contract, reader(doc.root), doc, decisionIds ?? new Set())];
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
