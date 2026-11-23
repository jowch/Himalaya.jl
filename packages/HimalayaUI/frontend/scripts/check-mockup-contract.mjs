#!/usr/bin/env node
import { readFileSync, existsSync } from "node:fs";
import { join, dirname } from "node:path";
import { fileURLToPath } from "node:url";

const FRONTEND = join(dirname(fileURLToPath(import.meta.url)), "..");
const CONTRACT = join(FRONTEND, "..", "..", "..", "docs", "mockup-contract.md");
const NAME = "docs/mockup-contract.md";
const STATUSES = new Set(["built", "untested", "superseded", "planned"]);
const ID = /^[a-z0-9]+(?:-[a-z0-9]+)*$/;

function cells(line) {
  return line.trim().replace(/^\|/, "").replace(/\|$/, "").split("|").map((c) => c.trim());
}

export function parseContract(text) {
  const lines = text.split("\n");
  const version = /Mock-up version: `([^`]+)`/.exec(text)?.[1];
  const sha256 = /Mock-up sha256: `([0-9a-f]{64})`/.exec(text)?.[1];
  const rows = [];
  lines.forEach((line, i) => {
    if (!line.trim().startsWith("|")) return;
    const [id, behaviour, status, evidence] = cells(line);
    if (!id || !ID.test(id) || id === "id") return;
    rows.push({ id, behaviour: behaviour ?? "", status: status ?? "", evidence: evidence ?? "", line: i + 1 });
  });
  return { version, sha256, rows };
}

function builtProblem(row, read) {
  const m = /^(\S+) › (.+)$/.exec(row.evidence);
  if (!m) return "a built row's evidence must be `path › test title`";
  const [, path, title] = m;
  const content = read(path);
  if (content === undefined) return `${path} does not exist`;
  if (!content.includes(title)) return `${path} does not contain the test "${title}"`;
  return null;
}

export function checkContract(contract, read) {
  const problems = [];
  if (!contract.version) problems.push(`${NAME} does not pin the mock-up version`);
  if (!contract.sha256) problems.push(`${NAME} does not pin the mock-up sha256`);
  const seen = new Set();
  for (const row of contract.rows) {
    const at = `${NAME}:${row.line} ${row.id}`;
    if (seen.has(row.id)) problems.push(`${at}: the id is used more than once`);
    seen.add(row.id);
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

function readFrontend(path) {
  const full = join(FRONTEND, path);
  return existsSync(full) ? readFileSync(full, "utf8") : undefined;
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  const contract = parseContract(readFileSync(CONTRACT, "utf8"));
  const problems = checkContract(contract, readFrontend);
  const counts = {};
  for (const r of contract.rows) counts[r.status] = (counts[r.status] ?? 0) + 1;
  if (problems.length) {
    for (const p of problems) console.error(p);
    console.error(`\n${problems.length} problem(s) in ${NAME}. Each row needs a status and the evidence that status requires; see the Statuses table at the top of the file.`);
    process.exit(1);
  }
  console.log(`${NAME}: ${contract.rows.length} behaviours (${Object.entries(counts).map(([k, v]) => `${v} ${k}`).join(", ")}) against mock-up ${contract.version}`);
}
