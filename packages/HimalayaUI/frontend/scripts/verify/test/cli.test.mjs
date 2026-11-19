import { test } from "node:test";
import assert from "node:assert/strict";
import { main, COMMANDS, helpText } from "../cli.mjs";
import { loadConfig } from "../config.mjs";

function sink() {
  const chunks = [];
  return { write: (s) => chunks.push(s), text: () => chunks.join(""), json: () => JSON.parse(chunks.join("")) };
}

const cfg = loadConfig({ DEMO_API_PORT: "1", DEMO_UI_PORT: "2", VERIFY_CDP_PORT: "3", VERIFY_STATE_DIR: "/tmp/verify-cli-test" });

test("help lists every command", async () => {
  const out = sink();
  assert.equal(await main(["help"], { cfg, out }), 0);
  for (const name of Object.keys(COMMANDS)) assert.match(out.text(), new RegExp(name));
});

test("per-command help prints usage", () => {
  assert.match(helpText("click"), /--role <role>/);
});

test("unknown command fails with a hint", async () => {
  const out = sink();
  assert.equal(await main(["frobnicate"], { cfg, out }), 1);
  assert.match(out.json().hint, /bin\/verify help/);
});

test("bad options fail as JSON with a hint", async () => {
  const out = sink();
  assert.equal(await main(["click", "--colour", "red"], { cfg, out }), 1);
  const r = out.json();
  assert.equal(r.ok, false);
  assert.match(r.error, /unknown option --colour/);
});

test("api --dry-run prints the request without sending it", async () => {
  const out = sink();
  assert.equal(await main(["api", "patch", "/api/exposures/3/status", "--data", '{"status":"rejected"}', "--dry-run"], { cfg, out }), 0);
  const r = out.json();
  assert.equal(r.dryRun, true);
  assert.equal(r.method, "PATCH");
  assert.deepEqual(r.body, { status: "rejected" });
});

test("api rejects bad JSON and non-/api paths", async () => {
  let out = sink();
  assert.equal(await main(["api", "POST", "/api/x", "--data", "{nope", "--dry-run"], { cfg, out }), 1);
  assert.match(out.json().error, /not valid JSON/);
  out = sink();
  assert.equal(await main(["api", "GET", "/health"], { cfg, out }), 1);
  assert.match(out.json().error, /must start with \/api\//);
});

test("api against a dead backend says how to fix it", async () => {
  const out = sink();
  assert.equal(await main(["api", "GET", "/api/health"], { cfg, out }), 1);
  assert.match(out.json().hint, /seed/);
});

test("seed --dry-run lists the steps, including deletes under --reset", async () => {
  const out = sink();
  assert.equal(await main(["seed", "--reset", "--dry-run"], { cfg, out }), 0);
  const steps = out.json().steps.join("\n");
  assert.match(steps, /delete .*demo/);
  assert.match(steps, /demo\.sh detached/);
});

test("doctor reports dead servers with fixes", async () => {
  const out = sink();
  assert.equal(await main(["doctor"], { cfg, out }), 1);
  const backend = out.json().checks.find((c) => c.name === "backend");
  assert.equal(backend.ok, false);
  assert.match(backend.fix, /bin\/verify seed/);
});

test("page commands without an open browser say to open one", async () => {
  const out = sink();
  assert.equal(await main(["snapshot"], { cfg, out }), 1);
  assert.match(out.json().hint, /bin\/verify open/);
});

test("missing positional argument fails with usage", async () => {
  const out = sink();
  assert.equal(await main(["open"], { cfg, out }), 1);
  assert.match(out.json().hint, /usage: bin\/verify open <route>/);
});
