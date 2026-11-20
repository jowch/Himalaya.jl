import { test } from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { main } from "../cli.mjs";
import { loadConfig } from "../config.mjs";
import { stubBackend } from "./stub.mjs";

function sink() {
  const chunks = [];
  return { write: (s) => chunks.push(s), json: () => JSON.parse(chunks.join("")) };
}

async function seedEmpty(experiments) {
  const stub = await stubBackend({ "/api/health": { status: "ok" }, "/api/experiments": experiments });
  const state = fs.mkdtempSync(path.join(os.tmpdir(), "verify-seed-"));
  const cfg = loadConfig({ ...stub.env, VERIFY_STATE_DIR: state, VERIFY_CDP_PORT: "1" });
  const out = sink();
  try {
    const code = await main(["seed", "--empty"], { cfg, out });
    return { code, r: out.json() };
  } finally {
    await stub.close();
  }
}

test("seed --empty refuses an app that already has experiments", async () => {
  const { code, r } = await seedEmpty([{ id: 1, data_dir: "/x/data", ingest_status: "complete" }]);
  assert.equal(code, 1);
  assert.match(r.error, /needs an empty DB.*1 experiment/);
  assert.match(r.hint, /DEMO_API_PORT=8081/);
});

test("seed --empty accepts an app with no experiments as already up", async () => {
  const { code, r } = await seedEmpty([]);
  assert.equal(code, 0, JSON.stringify(r));
  assert.equal(r.alreadyUp, true);
});

test("stop treats a pid file that no longer names demo.sh as stale and leaves the process alone", async () => {
  const state = fs.mkdtempSync(path.join(os.tmpdir(), "verify-stop-"));
  fs.writeFileSync(path.join(state, "demo.pid"), JSON.stringify({ pid: process.pid }));
  const cfg = loadConfig({ VERIFY_STATE_DIR: state });
  const out = sink();
  assert.equal(await main(["stop"], { cfg, out }), 0);
  const r = out.json();
  assert.equal(r.stopped, false);
  assert.match(r.reason, /not a demo\.sh process/);
  assert.equal(fs.existsSync(path.join(state, "demo.pid")), false);
});

test("--timeout must be an integer number of ms", async () => {
  const out = sink();
  const cfg = loadConfig({ VERIFY_STATE_DIR: os.tmpdir() });
  assert.equal(await main(["seed", "--timeout", "10m"], { cfg, out }), 1);
  assert.match(out.json().error, /--timeout must be a non-negative integer/);
});
