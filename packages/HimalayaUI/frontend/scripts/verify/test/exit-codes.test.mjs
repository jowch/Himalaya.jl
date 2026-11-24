import { test } from "node:test";
import assert from "node:assert/strict";
import { execFileSync } from "node:child_process";
import { readFileSync } from "node:fs";
import { join, dirname } from "node:path";
import { fileURLToPath } from "node:url";

const REPO = join(dirname(fileURLToPath(import.meta.url)), "..", "..", "..", "..", "..", "..");
const SHELL_OR = /\|\|\s*(?=\{|:(?:\s|;|$)|[A-Za-z_][\w-]*(?:[\s=;]|$))/g;
const FAILS = /\bexit\s+(?:[1-9]\d*\b|"?\$\?"?)|\breturn\s+[1-9]\d*\b|\bfailed=1\b|\bfailed\+=\(/;
const CLEANUP = /^\s*(?:trap\b|.*\bkill\b|.*\bwait\b|.*\bbin\/verify\s+(?:stop|close)\b|.*\bgit\s+worktree\s+remove\b)/;
const ERREXIT_OFF = /\bset\s+(?:[+][a-zA-Z]*e|\+o\s+(?:errexit|pipefail))\b/;
const MAKE_IGNORE = /^\t-/;

const tracked = execFileSync("git", ["ls-files", "-z"], { cwd: REPO, encoding: "utf8" }).split("\0").filter(Boolean);
const shellLike = tracked.filter((f) => f.endsWith(".sh") || f.startsWith("bin/") || f === "Makefile" || /^\.github\/workflows\/.*\.ya?ml$/.test(f));
const verifyScripts = shellLike.filter((f) => f.endsWith(".sh") || f === "bin/verify");

const lines = (file) => readFileSync(join(REPO, file), "utf8").split("\n").map((text, i) => ({ file, line: i + 1, text }));

function handlers(text) {
  return [...text.matchAll(SHELL_OR)].filter((m) => !/\]\]\s*$/.test(text.slice(0, m.index))).map((m) => text.slice(m.index + 2));
}

export function ignoresFailure(text) {
  if (ERREXIT_OFF.test(text)) return true;
  return handlers(text).some((h) => !FAILS.test(h));
}

export function swallowedCheckTruth(entries) {
  return entries.filter((l) => /check-truth/.test(l.text) && (ignoresFailure(l.text) || /check-truth[^|]*\|(?!\|)/.test(l.text) || MAKE_IGNORE.test(l.text)));
}

export function swallowedInVerifyScript(entries) {
  return entries.filter((l) => !CLEANUP.test(l.text) && ignoresFailure(l.text));
}

test("the guard accepts only a failure action after ||, and lets cleanup through", () => {
  const at = (text, file = "x.sh") => [{ file, line: 1, text }];
  const swallowed = [
    "bin/verify check-truth > f || true",
    "bin/verify check-truth || :",
    "bin/verify check-truth || exit 0",
    "bin/verify check-truth > f || status=$?",
    "bin/verify check-truth > f || rc=1",
    "bin/verify check-truth > f || echo broke",
    "bin/verify check-truth > f || { echo broke; }",
    'bin/verify check-truth > f || node "$help" broken f',
    "bin/verify check-truth | tee f",
    "set +e; bin/verify check-truth",
    "set +o errexit; bin/verify check-truth",
    "\t-bin/verify check-truth",
  ];
  for (const text of swallowed) assert.equal(swallowedCheckTruth(at(text)).length, 1, text);
  const failing = [
    'bin/verify check-truth --gate invariants > f || { node "$help" broken f; exit 1; }',
    'bin/verify check-truth > f || { echo "FAILED"; failed=1; }',
    "bin/verify check-truth > f || exit 1",
    'bin/verify check-truth > f || exit "$?"',
    "bin/verify check-truth > f",
  ];
  for (const text of failing) assert.equal(swallowedCheckTruth(at(text)).length, 0, text);
  assert.equal(swallowedInVerifyScript(at("trap 'bin/verify stop > /dev/null 2>&1 || true' EXIT")).length, 0);
  assert.equal(swallowedInVerifyScript(at('kill "$pid" 2>/dev/null || true')).length, 0);
  assert.equal(swallowedInVerifyScript(at("bin/verify snapshot > f || true")).length, 1);
  assert.equal(swallowedInVerifyScript(at('[[ -f $review ]] || review="$HISTORY"')).length, 0, "a failed test is a condition, not a swallowed exit");
  assert.equal(swallowedInVerifyScript(at("bin/verify snapshot > f || { sleep 1; continue; }")).length, 1);
  assert.equal(swallowedInVerifyScript(at('    if (!inv.checked || !inv.ok) { console.log(x); process.exit(1); }')).length, 0);
  assert.equal(swallowedInVerifyScript(at("    stale == 0 || error(\"stale\")")).length, 0);
});

test("no script or workflow ignores check-truth's exit code", () => {
  assert.ok(shellLike.includes("packages/HimalayaUI/scripts/verify_truth.sh"));
  assert.deepEqual(swallowedCheckTruth(shellLike.flatMap(lines)), []);
});

test("verify scripts ignore an exit code only to clean up", () => {
  assert.ok(verifyScripts.includes("bin/verify"));
  assert.deepEqual(swallowedInVerifyScript(verifyScripts.flatMap(lines)), []);
});

test("no workflow step or make recipe is allowed to fail", () => {
  const workflows = shellLike.filter((f) => f.startsWith(".github/workflows/"));
  assert.ok(workflows.length > 0);
  assert.deepEqual(workflows.flatMap(lines).filter((l) => /continue-on-error:\s*true/.test(l.text)), []);
  assert.deepEqual(lines("Makefile").filter((l) => MAKE_IGNORE.test(l.text)), []);
});
