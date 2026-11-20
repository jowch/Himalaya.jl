import { test } from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { main } from "../cli.mjs";
import { loadConfig } from "../config.mjs";
import { stubBackend } from "./stub.mjs";

const TRUTH = [
  "sample,bases,phases,lattice_A,auto,expected_peaks_q,images,note",
  "MO_H2O,MO_H2O_S1 MO_H2O_S2,Pn3m,102,Pn3m,,true,clean",
  'MO_PEG5,MO_PEG5_S3,Pn3m+Im3m,108+138,Im3m,,true,"coexistence, Pn3m hidden"',
].join("\n");

function routes(demoDir) {
  return {
    "/api/experiments": [
      { id: 3, data_dir: "/elsewhere/data", ingest_status: "complete" },
      { id: 4, data_dir: path.join(demoDir, "data"), ingest_status: "complete" },
    ],
    "/api/samples?experiment_id=4": [
      { id: 10, name: "MO_H2O_S1 (S01P01)" },
      { id: 11, name: "MO_PEG5_S3 (S01P02)" },
    ],
    "/api/experiments/4/export": [
      { name: "MO_H2O_S1 (S01P01)", exposures: [
        { filename: "MO_H2O_S1", selected: false, active_group_kind: "auto", indices: [{ phase: "Pn3m", basis: 0.087 }] },
        { filename: "MO_H2O_S2", selected: false, active_group_kind: "auto", indices: [{ phase: "Pn3m", basis: 0.087 }] },
      ] },
      { name: "MO_PEG5_S3 (S01P02)", exposures: [
        { filename: "MO_PEG5_S3", selected: false, active_group_kind: "auto", indices: [{ phase: "Im3m", basis: 0.064 }] },
      ] },
    ],
    "/api/samples/10/exposures": [
      { id: 2, filename: "MO_H2O_S2", status: null, selected: false, sample_id: 10 },
      { id: 1, filename: "MO_H2O_S1", status: null, selected: false, sample_id: 10 },
    ],
    "/api/samples/11/exposures": [{ id: 3, filename: "MO_PEG5_S3", status: null, selected: false, sample_id: 11 }],
    "/api/exposures/1/assignment": { exposure_id: 1, state: "indexed", members: [] },
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
    assert.deepEqual(h2o.auto_calls.map((c) => c.frame), ["MO_H2O_S1", "MO_H2O_S2"]);
    assert.equal(h2o.representative, "MO_H2O_S2");
    assert.equal(peg.assigned, "Im3m+Pn3m");
    assert.equal(peg.assigned_matches_truth, true);
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
