import { test } from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { main } from "../cli.mjs";
import { loadConfig } from "../config.mjs";
import { compareTraces, compareTruth } from "../truth.mjs";
import { heldTraces } from "../late.mjs";
import { stubBackend } from "./stub.mjs";

function sink() {
  const chunks = [];
  return { write: (s) => chunks.push(s), json: () => JSON.parse(chunks.join("")) };
}

const waiting = (filename) => ({ exposure_id: 1, filename, problem: "no_trace", reason: "no trace" });

function demoTree({ held = [] } = {}) {
  const demoDir = fs.mkdtempSync(path.join(os.tmpdir(), "verify-late-"));
  fs.mkdirSync(path.join(demoDir, "analysis", "automatic_analysis", "tot_files"), { recursive: true });
  if (held.length) {
    fs.mkdirSync(path.join(demoDir, "late"));
    for (const b of held) fs.writeFileSync(path.join(demoDir, "late", `${b}_tot.dat`), "0.1 1 0.1\n");
  }
  return demoDir;
}

test("compareTraces passes when exactly the held traces are reported as waiting", () => {
  assert.equal(compareTraces(["A_S1"], [waiting("A_S1")]).ok, true);
  assert.equal(compareTraces([], []).ok, true);
});

test("compareTraces fails on a delivered trace still waiting, a held one not reported, and a failed analysis", () => {
  const stale = compareTraces([], [waiting("A_S1")]);
  assert.equal(stale.ok, false);
  assert.deepEqual(stale.still_waiting, ["A_S1"]);
  const silent = compareTraces(["A_S1"], []);
  assert.equal(silent.ok, false);
  assert.deepEqual(silent.not_reported, ["A_S1"]);
  const failed = compareTraces([], [{ exposure_id: 2, filename: "B_S1", problem: "failed", reason: "bad header" }]);
  assert.equal(failed.ok, false);
  assert.deepEqual(failed.failed, [{ frame: "B_S1", reason: "bad header" }]);
});

test("a frame whose trace is held back must have no auto call yet", () => {
  const truth = [{ sample: "MO", bases: "MO_S1", phases: "Pn3m", auto: "Pn3m", drop: "", note: "" }];
  const app = (auto) => [{ id: 1, name: "MO_S1 (S01P01)", shown: { phase: null, state: "indexed", representative_id: 1 },
    frames: [{ id: 1, filename: "MO_S1", status: null, autoPhases: auto, assignment: { state: "indexed", phases: [] } }] }];
  assert.equal(compareTruth(truth, app([]), { held: new Set(["MO_S1"]) }).ok, true);
  const early = compareTruth(truth, app(["Pn3m"]), { held: new Set(["MO_S1"]) });
  assert.equal(early.ok, false);
  assert.match(early.samples[0].problems[0], /MO_S1=Pn3m \(expected none\)/);
  assert.equal(compareTruth(truth, app([])).ok, false);
});

test("heldTraces lists the bases waiting in late/", () => {
  assert.deepEqual(heldTraces(loadConfig({ DEMO_DIR: demoTree({ held: ["B_S2", "B_S1"] }) })), ["B_S1", "B_S2"]);
  assert.deepEqual(heldTraces(loadConfig({ DEMO_DIR: demoTree() })), []);
});

test("deliver-late moves the held traces into tot_files with a fresh mtime, then has nothing left", async () => {
  const demoDir = demoTree({ held: ["B_S1"] });
  const old = new Date("2026-01-01T00:00:00Z");
  fs.utimesSync(path.join(demoDir, "late", "B_S1_tot.dat"), old, old);
  const cfg = loadConfig({ DEMO_DIR: demoDir });
  const out = sink();
  assert.equal(await main(["deliver-late"], { cfg, out }), 0);
  assert.deepEqual(out.json().delivered, ["B_S1_tot.dat"]);
  const moved = path.join(demoDir, "analysis", "automatic_analysis", "tot_files", "B_S1_tot.dat");
  assert.ok(fs.statSync(moved).mtimeMs > old.getTime());
  assert.equal(fs.existsSync(path.join(demoDir, "late")), false);
  const again = sink();
  assert.equal(await main(["deliver-late"], { cfg, out: again }), 1);
  assert.match(again.json().hint, /seed --reset --late-trace/);
});

test("check-truth reads ingest problems and fails while a delivered trace is still waiting", async () => {
  const demoDir = demoTree({ held: ["MO_S1"] });
  fs.writeFileSync(path.join(demoDir, "TRUTH.csv"),
    "sample,bases,phases,lattice_A,auto,expected_peaks_q,images,drop,flare,note\nMO,MO_S1,Pn3m,100,Pn3m,,true,,,\n");
  const problems = { ingest_error: null, exposures: [waiting("MO_S1")] };
  const stub = await stubBackend({
    "/api/experiments": [{ id: 4, setups: [{ data_dir: path.join(demoDir, "data") }], ingest_status: "complete" }],
    "/api/samples?experiment_id=4": [{ id: 10, name: "MO_S1 (S01P01)", assignment_state: "indexed", phase: null, representative_exposure_id: 1 }],
    "/api/experiments/4/export": [{ name: "MO_S1 (S01P01)", exposures: [{ filename: "MO_S1", auto_phases: [] }] }],
    "/api/experiments/4/problems": () => problems,
    "/api/samples/10/exposures": [{ id: 1, filename: "MO_S1", status: null, sample_id: 10 }],
    "/api/exposures/1/assignment": { exposure_id: 1, state: "indexed", members: [] },
  });
  const cfg = loadConfig({ ...stub.env, DEMO_DIR: demoDir });
  try {
    const held = sink();
    assert.equal(await main(["check-truth"], { cfg, out: held }), 0, JSON.stringify(held.json(), null, 1));
    assert.deepEqual(held.json().traces.held, ["MO_S1"]);
    await main(["deliver-late"], { cfg, out: sink() });
    const delivered = sink();
    assert.equal(await main(["check-truth"], { cfg, out: delivered }), 1);
    assert.deepEqual(delivered.json().traces.still_waiting, ["MO_S1"]);
  } finally {
    await stub.close();
  }
});

test("seed --late-trace needs --reset, and the plan regenerates the demo with it", async () => {
  const cfg = loadConfig({ VERIFY_STATE_DIR: os.tmpdir() });
  const bare = sink();
  assert.equal(await main(["seed", "--late-trace"], { cfg, out: bare }), 1);
  assert.match(bare.json().error, /--late-trace needs --reset/);
  const plan = sink();
  assert.equal(await main(["seed", "--reset", "--late-trace", "--dry-run"], { cfg, out: plan }), 0);
  assert.ok(plan.json().steps.some((s) => s.includes("--late-trace")), JSON.stringify(plan.json().steps));
});

test("only a 404 on a held frame's trace is the app's expected answer while the trace is late", async () => {
  const { heldTrace404 } = await import("../views.mjs");
  const e = (url, status = 404, method = "GET") => ({ url: `http://127.0.0.1:5173${url}`, status, method });
  const held = new Set([11]);
  assert.equal(heldTrace404(e("/api/exposures/11/trace"), held), true);
  assert.equal(heldTrace404(e("/api/exposures/12/trace"), held), false);
  assert.equal(heldTrace404(e("/api/exposures/11/trace", 500), held), false);
  assert.equal(heldTrace404(e("/api/exposures/11/peaks"), held), false);
});

test("check-truth --ui tells the reader every held frame, since a neighbour's prefetch can ask for its trace", async () => {
  const { checkTruth } = await import("../truth.mjs");
  const demoDir = demoTree({ held: ["B_S2"] });
  fs.writeFileSync(path.join(demoDir, "TRUTH.csv"),
    "sample,bases,phases,lattice_A,auto,expected_peaks_q,images,drop,flare,note\nA,A_S1,none,,none,,true,,,\nB,B_S2,none,,none,,true,,,\n");
  const stub = await stubBackend({
    "/api/experiments": [{ id: 4, setups: [{ data_dir: path.join(demoDir, "data") }], ingest_status: "complete" }],
    "/api/samples?experiment_id=4": [
      { id: 10, name: "A_S1 (S01P01)", assignment_state: "indexed", phase: null, representative_exposure_id: 1 },
      { id: 11, name: "B_S2 (S01P02)", assignment_state: "indexed", phase: null, representative_exposure_id: 2 },
    ],
    "/api/experiments/4/export": [
      { name: "A_S1 (S01P01)", exposures: [{ filename: "A_S1", auto_phases: [] }] },
      { name: "B_S2 (S01P02)", exposures: [{ filename: "B_S2", auto_phases: [] }] },
    ],
    "/api/experiments/4/problems": { ingest_error: null, exposures: [waiting("B_S2")] },
    "/api/samples/10/exposures": [{ id: 1, filename: "A_S1", status: null, sample_id: 10 }],
    "/api/samples/11/exposures": [{ id: 2, filename: "B_S2", status: null, sample_id: 11 }],
    "/api/exposures/1/assignment": { exposure_id: 1, state: "indexed", members: [] },
    "/api/exposures/2/assignment": { exposure_id: 2, state: "indexed", members: [] },
  });
  const cfg = loadConfig({ ...stub.env, DEMO_DIR: demoDir });
  const asked = [];
  try {
    await checkTruth(cfg, { ui: true }, { views: async (_c, _e, samples) => { asked.push(...samples); return new Map(); } });
    assert.deepEqual(asked.map((s) => s.heldIds), [[2], [2]]);
    const read = (plot, notice) => async () => new Map([
      [10, { workspace: "A_S1", workspace_kept: "1/1", workspace_text: "Unreviewed", workspace_rail: { ticked: [], candidates: [], form_factor: false }, workspace_plot: null, workspace_notice: notice }],
      [11, { workspace: "B_S2", workspace_kept: "1/1", workspace_text: "Unreviewed", workspace_rail: { ticked: [], candidates: [], form_factor: false }, workspace_plot: plot, workspace_notice: notice }],
    ]);
    const good = await checkTruth(cfg, { ui: true }, { views: read("Waiting for this frame's trace.", "1 frame waiting for a trace") });
    assert.equal(good.invariants.plot_matches_trace.ok, true);
    assert.equal(good.traces.notice_ok, true);
    const stale = await checkTruth(cfg, { ui: true }, { views: read("No integrated trace for this frame.", null) });
    assert.equal(stale.ok, false);
    assert.deepEqual(stale.invariants.plot_matches_trace.violations, [{ sample: "B", frame: "B_S2", plot: "No integrated trace for this frame.", expected: "Waiting for this frame's trace." }]);
    assert.equal(stale.traces.notice_ok, false);
  } finally {
    await stub.close();
  }
});

test("plotMismatch: a held frame's plot says it is waiting, any other frame shows its trace", async () => {
  const { plotMismatch, WAITING_PLOT } = await import("../truth.mjs");
  const held = new Set(["B_S2"]);
  assert.equal(plotMismatch("B_S2", WAITING_PLOT, held), null);
  assert.deepEqual(plotMismatch("B_S2", "No integrated trace for this frame.", held), { frame: "B_S2", plot: "No integrated trace for this frame.", expected: WAITING_PLOT });
  assert.equal(plotMismatch("A_S1", null, held), null);
  assert.deepEqual(plotMismatch("A_S1", "No integrated trace for this frame.", held), { frame: "A_S1", plot: "No integrated trace for this frame.", expected: "the trace" });
  assert.equal(plotMismatch(null, null, held), null);
});

test("compareTraces checks the workspace notice when it was read", () => {
  assert.equal(compareTraces(["A_S1", "A_S2"], [waiting("A_S1"), waiting("A_S2")], "2 frames waiting for a trace").ok, true);
  assert.equal(compareTraces([], [], null).ok, true);
  const stale = compareTraces([], [], "1 frame waiting for a trace");
  assert.equal(stale.ok, false);
  assert.equal(stale.notice, "1 frame waiting for a trace");
  assert.equal(compareTraces(["A_S1"], [waiting("A_S1")], null).ok, false);
  assert.equal(compareTraces(["A_S1"], [waiting("A_S1")], undefined).ok, true);
});
