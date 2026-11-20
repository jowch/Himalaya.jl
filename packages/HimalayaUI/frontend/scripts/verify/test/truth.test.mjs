import { test } from "node:test";
import assert from "node:assert/strict";
import { compareTruth, normalizeCall, representativeFrame } from "../truth.mjs";

const TRUTH = [
  { sample: "MO_H2O", bases: "MO_H2O_S1 MO_H2O_S2", phases: "Pn3m", auto: "Pn3m", note: "" },
  { sample: "MO_PEG5", bases: "MO_PEG5_S3", phases: "Pn3m+Im3m", auto: "Im3m", note: "Pn3m √6 overlaps Im3m √10" },
];

function frame(filename, autoPhases, extra = {}) {
  return { id: Number(filename.at(-1)), filename, status: null, selected: false, autoPhases, assignment: { state: "indexed", phases: [] }, ...extra };
}

function app(overrides = {}) {
  return [
    { id: 1, name: "MO_H2O_S1 (S01P01)", frames: overrides.h2o ?? [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", ["Pn3m"])] },
    { id: 2, name: "MO_PEG5_S3 (S01P02)", frames: overrides.peg ?? [frame("MO_PEG5_S3", ["Im3m"])] },
  ];
}

test("normalizeCall sorts, dedupes and maps empty to none", () => {
  assert.equal(normalizeCall("Pn3m+Im3m"), "Im3m+Pn3m");
  assert.equal(normalizeCall(["Pn3m", "Pn3m"]), "Pn3m");
  assert.equal(normalizeCall([]), "none");
  assert.equal(normalizeCall("none"), "none");
});

test("matching auto calls pass and the truth gap is reported, not failed", () => {
  const r = compareTruth(TRUTH, app());
  assert.equal(r.ok, true);
  assert.deepEqual(r.summary.known_gaps, ["MO_PEG5"]);
  assert.equal(r.summary.frames_found, 3);
});

test("--strict fails where the auto call misses the generated phases", () => {
  const r = compareTruth(TRUTH, app(), { strict: true });
  assert.equal(r.ok, false);
  assert.deepEqual(r.summary.failed, ["MO_PEG5"]);
  assert.match(r.samples[1].problems[0], /misses the generated truth Im3m\+Pn3m/);
});

test("a changed auto call is a regression", () => {
  const r = compareTruth(TRUTH, app({ h2o: [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", ["Ia3d"])] }));
  assert.equal(r.ok, false);
  assert.match(r.samples[0].problems[0], /MO_H2O_S2=Ia3d, expected Pn3m/);
});

test("a missing frame and a missing sample fail", () => {
  const r = compareTruth(TRUTH, [app()[0]].map((s) => ({ ...s, frames: s.frames.slice(0, 1) })));
  assert.equal(r.ok, false);
  assert.deepEqual(r.samples[0].frames.missing, ["MO_H2O_S2"]);
  assert.deepEqual(r.samples[1].problems, ["sample missing from the app"]);
});

test("an app sample mixing two truth samples is a stray", () => {
  const merged = [{ id: 9, name: "merged", frames: [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_PEG5_S3", ["Im3m"])] }];
  const r = compareTruth(TRUTH, merged);
  assert.equal(r.ok, false);
  assert.deepEqual(r.strays[0].owners.sort(), ["MO_H2O", "MO_PEG5"]);
});

test("representative frame: selected if kept, else newest kept, else none", () => {
  const a = frame("a1", []);
  const b = frame("b2", [], { selected: true });
  const c = frame("c3", []);
  assert.equal(representativeFrame([a, b, c]).filename, "b2");
  assert.equal(representativeFrame([a, { ...b, status: "rejected" }, c]).filename, "c3");
  assert.equal(representativeFrame([a, { ...c, status: "rejected" }]).filename, "a1");
  assert.equal(representativeFrame([{ ...a, status: "rejected" }]), null);
});

test("the assigned call is read from the representative frame and flagged when it sits elsewhere", () => {
  const assigned = { state: "indexed", phases: ["Pn3m", "Im3m"] };
  const onRep = compareTruth(TRUTH, app({ peg: [frame("MO_PEG5_S3", ["Im3m"], { assignment: assigned })] }));
  assert.equal(onRep.samples[1].assigned, "Im3m+Pn3m");
  assert.equal(onRep.samples[1].assigned_matches_truth, true);

  const offRep = compareTruth(TRUTH, app({ h2o: [frame("MO_H2O_S1", ["Pn3m"], { assignment: { state: "indexed", phases: ["Pn3m"] } }), frame("MO_H2O_S2", ["Pn3m"])] }));
  assert.equal(offRep.samples[0].assigned, null);
  assert.deepEqual(offRep.summary.assigned_off_representative, ["MO_H2O"]);
});

test("a form-factor assignment reads as none", () => {
  const ff = { state: "form_factor", phases: [] };
  const r = compareTruth(TRUTH, app({ peg: [frame("MO_PEG5_S3", ["Im3m"], { assignment: ff })] }));
  assert.equal(r.samples[1].assigned, "none");
  assert.equal(r.samples[1].assigned_matches_truth, false);
});
