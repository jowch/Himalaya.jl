import { test } from "node:test";
import assert from "node:assert/strict";
import { compareTruth, normalizeCall, shownStatus } from "../truth.mjs";

const TRUTH = [
  { sample: "MO_H2O", bases: "MO_H2O_S1 MO_H2O_S2", phases: "Pn3m", auto: "Pn3m", drop: "", note: "" },
  { sample: "MO_PEG5", bases: "MO_PEG5_S3", phases: "Pn3m+Im3m", auto: "Im3m", drop: "", note: "Pn3m √6 overlaps Im3m √10" },
];

const NOT_INDEXED = { phase: null, state: "indexed" };

function frame(filename, autoPhases, extra = {}) {
  return { id: Number(filename.at(-1)), filename, status: null, selected: false, autoPhases, assignment: { state: "indexed", phases: [] }, ...extra };
}

function sample(id, name, frames, shown = {}) {
  return { id, name, frames, shown: { ...NOT_INDEXED, representative_id: frames.at(-1).id, ...shown } };
}

function app(overrides = {}) {
  return [
    sample(1, "MO_H2O_S1 (S01P01)", overrides.h2o ?? [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", ["Pn3m"])], overrides.h2oShown),
    sample(2, "MO_PEG5_S3 (S01P02)", overrides.peg ?? [frame("MO_PEG5_S3", ["Im3m"])], overrides.pegShown),
  ];
}

const dropped = { status: "rejected" };

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
  assert.match(r.samples[0].problems[0], /MO_H2O_S2=Ia3d \(expected Pn3m\)/);
});

test("a missing frame and a missing sample fail", () => {
  const r = compareTruth(TRUTH, [app()[0]].map((s) => ({ ...s, frames: s.frames.slice(0, 1) })));
  assert.equal(r.ok, false);
  assert.deepEqual(r.samples[0].frames.missing, ["MO_H2O_S2"]);
  assert.deepEqual(r.samples[1].problems, ["sample missing from the app"]);
});

test("an app sample mixing two truth samples is a stray", () => {
  const merged = [sample(9, "merged", [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_PEG5_S3", ["Im3m"])])];
  const r = compareTruth(TRUTH, merged);
  assert.equal(r.ok, false);
  assert.deepEqual(r.strays[0].owners.sort(), ["MO_H2O", "MO_PEG5"]);
});

test("a frame in TRUTH's drop column must call none, and whether it is dropped is reported", () => {
  const truth = [{ ...TRUTH[0], drop: "MO_H2O_S2" }, TRUTH[1]];
  const kept = compareTruth(truth, app({ h2o: [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", [])] }));
  assert.equal(kept.ok, true);
  assert.deepEqual(kept.samples[0].frames.bad_still_kept, ["MO_H2O_S2"]);
  assert.deepEqual(kept.summary.bad_still_kept, ["MO_H2O_S2"]);

  const triaged = compareTruth(truth, app({ h2o: [frame("MO_H2O_S1", ["Pn3m"], dropped), frame("MO_H2O_S2", [], dropped)] }));
  assert.deepEqual(triaged.samples[0].frames.bad_still_kept, []);
  assert.deepEqual(triaged.samples[0].frames.good_dropped, ["MO_H2O_S1"]);

  const indexed = compareTruth(truth, app({ h2o: [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", ["Pn3m"])] }));
  assert.equal(indexed.ok, false);
  assert.match(indexed.samples[0].problems[0], /MO_H2O_S2=Pn3m \(expected none\)/);
});

test("the shown status is what GET /api/samples gives the corpus, on the frame it names", () => {
  assert.deepEqual(shownStatus({ phase: "Pn3m", state: "indexed" }), { status: "phase", call: "Pn3m" });
  assert.deepEqual(shownStatus({ phase: null, state: "form_factor" }), { status: "form_factor", call: "none" });
  assert.deepEqual(shownStatus(NOT_INDEXED), { status: "unreviewed", call: null });

  const r = compareTruth(TRUTH, app({ pegShown: { phase: "Im3m" } }));
  assert.deepEqual(r.samples[1].shown, { status: "phase", call: "Im3m", frame: "MO_PEG5_S3" });
  assert.equal(r.samples[1].shown_matches_truth, false);
  assert.equal(r.samples[0].shown_matches_truth, null);
  assert.equal(r.summary.shown_reviewed, 1);
});

test("a broken invariant fails the gate; an unchecked one does not", () => {
  const r = compareTruth(TRUTH, app());
  assert.equal(r.ok, true);
  assert.deepEqual(r.summary.invariants_failed, []);
  const broken = compareTruth(TRUTH, app({ pegShown: { phase: "Im3m" } }));
  assert.equal(broken.ok, false);
  assert.deepEqual(broken.summary.invariants_failed, ["status_matches_assignment"]);
  assert.equal(r.invariants.views_agree.checked, false);
  assert.equal(r.invariants.views_agree.ok, null);
  assert.equal("invariants" in r.samples[0], false);
});

test("status_matches_assignment: unconfirmed reads unreviewed, confirmed reads that call", () => {
  const assigned = { state: "indexed", phases: ["Pn3m"] };
  const h2oAssigned = [frame("MO_H2O_S1", ["Pn3m"], { assignment: assigned }), frame("MO_H2O_S2", ["Pn3m"])];

  const fresh = compareTruth(TRUTH, app());
  assert.equal(fresh.invariants.status_matches_assignment.ok, true, "an auto call alone is not a confirmed call");

  const hidden = compareTruth(TRUTH, app({ h2o: h2oAssigned }));
  assert.deepEqual(hidden.invariants.status_matches_assignment.violations, [
    { sample: "MO_H2O", shown: "unreviewed", confirmed: [{ frame: "MO_H2O_S1", call: "Pn3m" }] },
  ]);

  const shown = compareTruth(TRUTH, app({ h2o: h2oAssigned, h2oShown: { phase: "Pn3m" } }));
  assert.equal(shown.invariants.status_matches_assignment.ok, true);

  const wrong = compareTruth(TRUTH, app({ h2o: h2oAssigned, h2oShown: { phase: "Ia3d" } }));
  assert.equal(wrong.invariants.status_matches_assignment.ok, false);

  const unconfirmedShown = compareTruth(TRUTH, app({ pegShown: { phase: "Im3m" } }));
  assert.deepEqual(unconfirmedShown.invariants.status_matches_assignment.violations, [{ sample: "MO_PEG5", shown: "Im3m", confirmed: [] }]);

  const coexistence = { state: "indexed", phases: ["Pn3m", "Im3m"] };
  const peg = compareTruth(TRUTH, app({ peg: [frame("MO_PEG5_S3", ["Im3m"], { assignment: coexistence })], pegShown: { phase: "Pn3m" } }));
  assert.equal(peg.invariants.status_matches_assignment.ok, true, "the corpus shows one phase of a coexistence call");

  const onDropped = compareTruth(TRUTH, app({ h2o: [frame("MO_H2O_S1", ["Pn3m"], { assignment: assigned, ...dropped }), frame("MO_H2O_S2", ["Pn3m"])] }));
  assert.equal(onDropped.invariants.status_matches_assignment.ok, true);

  const ff = { state: "form_factor", phases: [] };
  const ffHidden = compareTruth(TRUTH, app({ peg: [frame("MO_PEG5_S3", ["Im3m"], { assignment: ff })] }));
  assert.equal(ffHidden.invariants.status_matches_assignment.ok, false);
  const ffShown = compareTruth(TRUTH, app({ peg: [frame("MO_PEG5_S3", ["Im3m"], { assignment: ff })], pegShown: { state: "form_factor" } }));
  assert.equal(ffShown.invariants.status_matches_assignment.ok, true);
});

test("representative_kept: the corpus status describes a dropped frame while a kept one exists", () => {
  const h2o = [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", ["Pn3m"], dropped)];
  const r = compareTruth(TRUTH, app({ h2o }));
  assert.deepEqual(r.invariants.representative_kept.violations, [
    { sample: "MO_H2O", view: "corpus_api", frame: "MO_H2O_S2", kept: ["MO_H2O_S1"] },
  ]);

  const allDropped = compareTruth(TRUTH, app({ h2o: h2o.map((f) => ({ ...f, ...dropped })) }));
  assert.equal(allDropped.invariants.representative_kept.ok, true);

  const onKept = compareTruth(TRUTH, app({ h2o, h2oShown: { representative_id: 1 } }));
  assert.equal(onKept.invariants.representative_kept.ok, true);
});

function withViews(samples, views, corpusText = "Not indexed") {
  return samples.map((s) => ({ ...s, views: views[s.id] ?? { loupe: s.frames.at(-1).filename, focus: s.frames.at(-1).filename }, corpusText }));
}

test("views_agree compares the corpus frame with the frames Loupe and Focus open on", () => {
  const agree = compareTruth(TRUTH, withViews(app(), {}));
  assert.equal(agree.invariants.views_agree.checked, true);
  assert.equal(agree.invariants.views_agree.ok, true);

  const split = compareTruth(TRUTH, withViews(app(), { 1: { loupe: "MO_H2O_S1", focus: "MO_H2O_S1" } }));
  assert.deepEqual(split.invariants.views_agree.violations, [
    { sample: "MO_H2O", corpus_api: "MO_H2O_S2", loupe: "MO_H2O_S1", focus: "MO_H2O_S1" },
  ]);
  assert.deepEqual(split.samples[0].views, { loupe: "MO_H2O_S1", focus: "MO_H2O_S1" });
});

test("a view opening on a dropped frame breaks representative_kept", () => {
  const h2o = [frame("MO_H2O_S1", ["Pn3m"], dropped), frame("MO_H2O_S2", ["Pn3m"])];
  const r = compareTruth(TRUTH, withViews(app({ h2o }), { 1: { loupe: "MO_H2O_S1", focus: "MO_H2O_S2" } }));
  assert.deepEqual(r.invariants.representative_kept.violations, [
    { sample: "MO_H2O", view: "loupe", frame: "MO_H2O_S1", kept: ["MO_H2O_S2"] },
  ]);
});

test("with --ui, the corpus text must match the status read from the API", () => {
  const same = compareTruth(TRUTH, withViews(app(), {}));
  assert.equal(same.ok, true);
  assert.equal(same.samples[0].shown.corpus_text, "Not indexed");

  const differ = compareTruth(TRUTH, withViews(app({ pegShown: { phase: "Im3m" } }), {}));
  assert.equal(differ.ok, false);
  assert.match(differ.samples[1].problems[0], /corpus reads "Not indexed" but GET \/api\/samples says phase/);
});

test("a Loupe or Focus read that finds no frame is a problem, and the view invariants report unchecked", () => {
  const r = compareTruth(TRUTH, withViews(app(), { 1: { loupe: null, focus: "MO_H2O_S2" }, 2: { loupe: "MO_PEG5_S3", focus: null } }));
  assert.equal(r.ok, false);
  assert.match(r.samples[0].problems[0], /could not read which frame Loupe opens on/);
  assert.match(r.samples[1].problems[0], /could not read which frame Focus opens on/);
  for (const name of ["representative_kept", "views_agree"]) {
    assert.equal(r.invariants[name].checked, false);
    assert.equal(r.invariants[name].ok, null);
  }
  assert.deepEqual(r.invariants.views_agree.violations, []);
  assert.equal(r.invariants.status_matches_assignment.checked, true);
});

test("Focus opening nothing is expected when every frame is dropped", () => {
  const h2o = [frame("MO_H2O_S1", ["Pn3m"], dropped), frame("MO_H2O_S2", ["Pn3m"], dropped)];
  const r = compareTruth(TRUTH, withViews(app({ h2o }), { 1: { loupe: "MO_H2O_S1", focus: null } }));
  assert.deepEqual(r.samples[0].problems, []);
  assert.equal(r.invariants.views_agree.checked, true);
});
