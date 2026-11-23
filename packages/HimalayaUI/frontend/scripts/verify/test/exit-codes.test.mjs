import { test } from "node:test";
import assert from "node:assert/strict";
import { execFileSync } from "node:child_process";
import { readFileSync } from "node:fs";
import { join, dirname } from "node:path";
import { fileURLToPath } from "node:url";

const REPO = join(dirname(fileURLToPath(import.meta.url)), "..", "..", "..", "..", "..", "..");
const SWALLOW = /\|\|\s*(?:true\b|:(?=\s|;|$)|exit\s+0\b)|;\s*true\s*$|\|\|\s*\w+=\$\?|\bset\s+\+e\b/;
const CLEANUP = /^\s*(?:trap\b|.*\bkill\b|.*\bwait\b|.*\bbin\/verify\s+(?:stop|close)\b|.*\bgit\s+worktree\s+remove\b)/;

const tracked = execFileSync("git", ["ls-files", "-z"], { cwd: REPO, encoding: "utf8" }).split("\0").filter(Boolean);
const shellLike = tracked.filter((f) => f.endsWith(".sh") || f.startsWith("bin/") || f === "Makefile" || /^\.github\/workflows\/.*\.ya?ml$/.test(f));
const verifyScripts = shellLike.filter((f) => f.endsWith(".sh") || f === "bin/verify");

const lines = (file) => readFileSync(join(REPO, file), "utf8").split("\n").map((text, i) => ({ file, line: i + 1, text }));

export function swallowedCheckTruth(entries) {
  return entries.filter((l) => /check-truth/.test(l.text) && SWALLOW.test(l.text));
}

export function swallowedInVerifyScript(entries) {
  return entries.filter((l) => SWALLOW.test(l.text) && !CLEANUP.test(l.text));
}

test("the guard spots each way of ignoring an exit code, and lets cleanup through", () => {
  const at = (text) => [{ file: "x.sh", line: 1, text }];
  for (const text of ["bin/verify check-truth > f || true", "bin/verify check-truth || :", "bin/verify check-truth || exit 0", "bin/verify check-truth > f || status=$?", "set +e; bin/verify check-truth"]) {
    assert.equal(swallowedCheckTruth(at(text)).length, 1, text);
  }
  assert.equal(swallowedCheckTruth(at('bin/verify check-truth --gate invariants > f || node "$help" broken f')).length, 0);
  assert.equal(swallowedInVerifyScript(at("trap 'bin/verify stop > /dev/null 2>&1 || true' EXIT")).length, 0);
  assert.equal(swallowedInVerifyScript(at('kill "$pid" 2>/dev/null || true')).length, 0);
  assert.equal(swallowedInVerifyScript(at("bin/verify snapshot > f || true")).length, 1);
});

test("no script or workflow ignores check-truth's exit code", () => {
  assert.ok(shellLike.includes("packages/HimalayaUI/scripts/verify_truth.sh"));
  assert.deepEqual(swallowedCheckTruth(shellLike.flatMap(lines)), []);
});

test("verify scripts ignore an exit code only to clean up", () => {
  assert.ok(verifyScripts.includes("bin/verify"));
  assert.deepEqual(swallowedInVerifyScript(verifyScripts.flatMap(lines)), []);
});

test("no workflow step is allowed to fail", () => {
  const workflows = shellLike.filter((f) => f.startsWith(".github/workflows/"));
  assert.ok(workflows.length > 0);
  assert.deepEqual(workflows.flatMap(lines).filter((l) => /continue-on-error:\s*true/.test(l.text)), []);
});
