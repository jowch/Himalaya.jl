import { test } from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { main } from "../cli.mjs";
import { checkTruth } from "../truth.mjs";
import { loadConfig } from "../config.mjs";
import { stubBackend } from "./stub.mjs";

const TRUTH = [
  "sample,bases,phases,lattice_A,auto,expected_peaks_q,images,drop,flare,note",
  "MO_H2O,MO_H2O_S1 MO_H2O_S2,Pn3m,102,Pn3m,,true,MO_H2O_S2,,clean",
  'MO_PEG5,MO_PEG5_S3,Pn3m+Im3m,108+138,Im3m,,true,,,"coexistence, Pn3m hidden"',
].join("\n");

function routes(demoDir) {
  return {
    "/api/experiments": [
      { id: 3, data_dir: "/elsewhere/data", ingest_status: "complete" },
      { id: 4, data_dir: path.join(demoDir, "data"), ingest_status: "complete" },
    ],
    "/api/samples?experiment_id=4": [
      { id: 10, name: "MO_H2O_S1 (S01P01)", assignment_state: "indexed", phase: "Pn3m", representative_exposure_id: 1 },
      { id: 11, name: "MO_PEG5_S3 (S01P02)", assignment_state: "indexed", phase: "Pn3m", representative_exposure_id: 3 },
    ],
    "/api/experiments/4/export": [
      { name: "MO_H2O_S1 (S01P01)", exposures: [
        { filename: "MO_H2O_S1", selected: false, active_group_kind: "auto", indices: [{ phase: "Pn3m", basis: 0.087 }] },
        { filename: "MO_H2O_S2", selected: false, active_group_kind: "auto", indices: [] },
      ] },
      { name: "MO_PEG5_S3 (S01P02)", exposures: [
        { filename: "MO_PEG5_S3", selected: false, active_group_kind: "auto", indices: [{ phase: "Im3m", basis: 0.064 }] },
      ] },
    ],
    "/api/samples/10/exposures": [
      { id: 2, filename: "MO_H2O_S2", status: "rejected", selected: false, sample_id: 10 },
      { id: 1, filename: "MO_H2O_S1", status: null, selected: false, sample_id: 10 },
    ],
    "/api/samples/11/exposures": [{ id: 3, filename: "MO_PEG5_S3", status: null, selected: false, sample_id: 11 }],
    "/api/exposures/1/assignment": { exposure_id: 1, state: "indexed", members: [5] },
    "/api/exposures/1/indices": [{ id: 5, phase: "Pn3m", kind: "auto", status: "candidate" }],
    "/api/exposures/2/assignment": { exposure_id: 2, state: "indexed", members: [] },
    "/api/exposures/3/assignment": { exposure_id: 3, state: "indexed", members: [7, 9] },
    "/api/exposures/3/indices": [
      { id: 7, phase: "Pn3m", kind: "auto", status: "candidate" },
      { id: 8, phase: "Ia3d", kind: "auto", status: "candidate" },
      { id: 9, phase: "Im3m", kind: "auto", status: "candidate" },
    ],
  };
}

function sink() {
  const chunks = [];
  return { write: (s) => chunks.push(s), json: () => JSON.parse(chunks.join("")) };
}

test("check-truth reads the backend's real endpoint shapes end to end", async () => {
  const demoDir = fs.mkdtempSync(path.join(os.tmpdir(), "verify-truth-"));
  fs.writeFileSync(path.join(demoDir, "TRUTH.csv"), TRUTH);
  const stub = await stubBackend(routes(demoDir));
  const cfg = loadConfig({ ...stub.env, DEMO_DIR: demoDir });
  const out = sink();
  try {
    const code = await main(["check-truth"], { cfg, out });
    const r = out.json();
    assert.equal(code, 0, JSON.stringify(r, null, 1));
    assert.equal(r.experiment_id, 4);
    assert.equal(r.summary.frames_found, 3);
    assert.deepEqual(r.summary.known_gaps, ["MO_PEG5"]);
    const [h2o, peg] = r.samples;
    assert.deepEqual(h2o.auto_calls.map((c) => [c.frame, c.expected]), [["MO_H2O_S1", "Pn3m"], ["MO_H2O_S2", "none"]]);
    assert.deepEqual(h2o.shown, { status: "phase", call: "Pn3m", frame: "MO_H2O_S1" });
    assert.deepEqual(h2o.assigned_frames, [{ frame: "MO_H2O_S1", call: "Pn3m" }]);
    assert.deepEqual(peg.shown, { status: "phase", call: "Pn3m", frame: "MO_PEG5_S3" });
    assert.deepEqual(peg.assigned_frames, [{ frame: "MO_PEG5_S3", call: "Im3m+Pn3m" }]);
    assert.deepEqual(r.summary.invariants_failed, []);
    assert.ok(stub.seen.includes("/api/samples?experiment_id=4"));
    assert.ok(!stub.seen.some((u) => u.includes("experiment_id=3")));
  } finally {
    await stub.close();
  }
});

test("check-truth fails loudly when an endpoint stops answering", async () => {
  const demoDir = fs.mkdtempSync(path.join(os.tmpdir(), "verify-truth-"));
  fs.writeFileSync(path.join(demoDir, "TRUTH.csv"), TRUTH);
  const r0 = routes(demoDir);
  delete r0["/api/experiments/4/export"];
  const stub = await stubBackend(r0);
  const cfg = loadConfig({ ...stub.env, DEMO_DIR: demoDir });
  const out = sink();
  try {
    assert.equal(await main(["check-truth"], { cfg, out }), 1);
    assert.match(out.json().error, /export returned 404/);
  } finally {
    await stub.close();
  }
});

test("check-truth says to seed when the demo experiment is missing", async () => {
  const demoDir = fs.mkdtempSync(path.join(os.tmpdir(), "verify-truth-"));
  fs.writeFileSync(path.join(demoDir, "TRUTH.csv"), TRUTH);
  const stub = await stubBackend({ "/api/experiments": [] });
  const cfg = loadConfig({ ...stub.env, DEMO_DIR: demoDir });
  const out = sink();
  try {
    assert.equal(await main(["check-truth"], { cfg, out }), 1);
    assert.match(out.json().hint, /bin\/verify seed/);
  } finally {
    await stub.close();
  }
});

test("triage drops TRUTH's bad frames once and --undo keeps them again", async () => {
  const demoDir = fs.mkdtempSync(path.join(os.tmpdir(), "verify-truth-"));
  fs.writeFileSync(path.join(demoDir, "TRUTH.csv"), TRUTH);
  const r0 = routes(demoDir);
  r0["/api/samples/10/exposures"][0].status = null;
  r0["PATCH /api/exposures/2/status"] = { id: 2, status: "rejected" };
  const stub = await stubBackend(r0);
  const cfg = loadConfig({ ...stub.env, DEMO_DIR: demoDir });
  try {
    const dry = sink();
    assert.equal(await main(["triage", "--dry-run"], { cfg, out: dry }), 0);
    assert.deepEqual(dry.json().changed.map((c) => c.frame), ["MO_H2O_S2"]);
    assert.deepEqual(stub.sent, []);

    const out = sink();
    assert.equal(await main(["triage"], { cfg, out }), 0);
    assert.deepEqual(stub.sent, [{ key: "PATCH /api/exposures/2/status", body: { status: "rejected" } }]);

    r0["/api/samples/10/exposures"][0].status = "rejected";
    const again = sink();
    assert.equal(await main(["triage"], { cfg, out: again }), 0);
    assert.deepEqual(again.json().changed, []);

    const undo = sink();
    assert.equal(await main(["triage", "--undo"], { cfg, out: undo }), 0);
    assert.deepEqual(stub.sent.at(-1), { key: "PATCH /api/exposures/2/status", body: { status: null } });
  } finally {
    await stub.close();
  }
});

test("check-truth --ui hands the reader each sample's exposure order and folds in what it saw", async () => {
  const demoDir = fs.mkdtempSync(path.join(os.tmpdir(), "verify-truth-"));
  fs.writeFileSync(path.join(demoDir, "TRUTH.csv"), TRUTH);
  const stub = await stubBackend(routes(demoDir));
  const cfg = loadConfig({ ...stub.env, DEMO_DIR: demoDir });
  const asked = [];
  const views = async (_cfg, experimentId, samples) => {
    asked.push({ experimentId, samples });
    return new Map([
      [10, { workspace_kept: "2/2", workspace: "MO_H2O_S2", workspace_text: "Pn3m" }],
      [11, { workspace_kept: "1/1", workspace: "MO_PEG5_S3", workspace_text: "Pn3m" }],
    ]);
  };
  try {
    const r = await checkTruth(cfg, { ui: true }, { views });
    assert.equal(asked[0].experimentId, 4);
    assert.deepEqual(asked[0].samples[0], { id: 10, name: "MO_H2O_S1 (S01P01)", frameOrder: ["MO_H2O_S2", "MO_H2O_S1"] });
    assert.equal(r.ok, false);
    assert.deepEqual(r.summary.invariants_failed, ["representative_kept", "views_agree"]);
    assert.equal(r.invariants.views_agree.checked, true);
    assert.deepEqual(r.invariants.views_agree.violations, [{ sample: "MO_H2O", corpus_api: "MO_H2O_S1", workspace: "MO_H2O_S2" }]);
    assert.deepEqual(r.invariants.representative_kept.violations.map((v) => v.view), ["workspace"]);
  } finally {
    await stub.close();
  }
});

test("check-truth --ui fails when the workspace list, read in place after a write, lags the API", async () => {
  const demoDir = fs.mkdtempSync(path.join(os.tmpdir(), "verify-truth-"));
  fs.writeFileSync(path.join(demoDir, "TRUTH.csv"), TRUTH);
  const stub = await stubBackend(routes(demoDir));
  const cfg = loadConfig({ ...stub.env, DEMO_DIR: demoDir });
  const views = async () =>
    new Map([
      [10, { workspace: "MO_H2O_S1", workspace_text: "Unreviewed", workspace_reused: true }],
      [11, { workspace: "MO_PEG5_S3", workspace_text: "Pn3m", workspace_reused: true }],
    ]);
  try {
    const r = await checkTruth(cfg, { ui: true }, { views });
    assert.equal(r.workspace_read, "in place");
    assert.equal(r.ok, false);
    assert.ok(r.summary.invariants_failed.includes("workspace_status_matches"));
    assert.deepEqual(r.invariants.workspace_status_matches.violations, [{ sample: "MO_H2O", shown: "Pn3m", workspace: "Unreviewed" }]);
    assert.equal("workspace_reused" in r.samples[0].views, false);
  } finally {
    await stub.close();
  }
});
