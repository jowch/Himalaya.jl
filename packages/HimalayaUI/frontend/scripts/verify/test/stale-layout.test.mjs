import { test } from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { main } from "../cli.mjs";
import { loadConfig } from "../config.mjs";
import { staleDemoLayout } from "../demo.mjs";
import { stubBackend } from "./stub.mjs";

function sink() {
  const chunks = [];
  return { write: (s) => chunks.push(s), json: () => JSON.parse(chunks.join("")) };
}

const TRUTH = "sample,bases,phases,lattice_A,auto,expected_peaks_q,images,drop,flare,note\nMO,MO_S1,Pn3m,100,Pn3m,,true,,,\n";

function demo({ layout }) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), "verify-stale-"));
  fs.writeFileSync(path.join(dir, "TRUTH.csv"), TRUTH);
  fs.mkdirSync(path.join(dir, layout === "old" ? "data" : path.join("capillary", "data")), { recursive: true });
  return dir;
}

test("a demo with top-level data/ and no capillary/ is the layout from before the second setup", () => {
  assert.equal(staleDemoLayout(loadConfig({ DEMO_DIR: demo({ layout: "old" }) })), true);
  assert.equal(staleDemoLayout(loadConfig({ DEMO_DIR: demo({ layout: "new" }) })), false);
  assert.equal(staleDemoLayout(loadConfig({ DEMO_DIR: fs.mkdtempSync(path.join(os.tmpdir(), "verify-stale-")) })), false);
});

test("doctor fails demo-data on the old layout and points at seed --reset", async () => {
  const cfg = loadConfig({ DEMO_DIR: demo({ layout: "old" }), DEMO_API_PORT: "1", DEMO_UI_PORT: "1", VERIFY_CDP_PORT: "1" });
  const out = sink();
  assert.equal(await main(["doctor"], { cfg, out }), 1);
  const check = out.json().checks.find((c) => c.name === "demo-data");
  assert.equal(check.ok, false);
  assert.match(check.fix, /seed --reset/);
});

test("check-truth on the old layout says to reset instead of to seed", async () => {
  const stub = await stubBackend({ "/api/experiments": [] });
  const cfg = loadConfig({ ...stub.env, DEMO_DIR: demo({ layout: "old" }) });
  const out = sink();
  try {
    assert.equal(await main(["check-truth"], { cfg, out }), 1);
    assert.match(out.json().hint, /seed --reset/);
  } finally {
    await stub.close();
  }
});

test("seed refuses the old layout before starting anything, and --reset is the way out", async () => {
  const state = fs.mkdtempSync(path.join(os.tmpdir(), "verify-seed-"));
  const cfg = loadConfig({ DEMO_DIR: demo({ layout: "old" }), DEMO_API_PORT: "1", DEMO_UI_PORT: "1", VERIFY_CDP_PORT: "1", VERIFY_STATE_DIR: state });
  const out = sink();
  assert.equal(await main(["seed"], { cfg, out }), 1);
  assert.match(out.json().hint, /seed --reset/);
  const plan = sink();
  assert.equal(await main(["seed", "--reset", "--dry-run"], { cfg, out: plan }), 0);
});
