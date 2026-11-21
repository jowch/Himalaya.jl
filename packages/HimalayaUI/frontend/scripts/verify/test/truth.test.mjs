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

test("a flared frame is to be dropped, and its auto call may lose phases but not add one", () => {
  const truth = [{ ...TRUTH[0], drop: "MO_H2O_S2", flare: "MO_H2O_S2" }, TRUTH[1]];
  for (const call of [["Pn3m"], []]) {
    const r = compareTruth(truth, app({ h2o: [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", call)] }));
    assert.equal(r.ok, true);
    assert.deepEqual(r.samples[0].frames.bad_still_kept, ["MO_H2O_S2"]);
    assert.equal(r.samples[0].auto_calls[1].expected, "within Pn3m");
  }
  const invented = compareTruth(truth, app({ h2o: [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", ["Ia3d", "Pn3m"])] }));
  assert.equal(invented.ok, false);
  assert.match(invented.samples[0].problems[0], /MO_H2O_S2=Ia3d\+Pn3m \(expected within Pn3m\)/);
  const unflared = compareTruth(truth, app({ h2o: [frame("MO_H2O_S1", []), frame("MO_H2O_S2", ["Pn3m"])] }));
  assert.equal(unflared.ok, false);
});

test("--expect fresh fails on a frame already dropped; --expect triaged fails on a bad frame still kept", () => {
  const truth = [{ ...TRUTH[0], drop: "MO_H2O_S2" }, TRUTH[1]];
  const untouched = app({ h2o: [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", [])] });
  const triaged = app({ h2o: [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", [], dropped)], h2oShown: { representative_id: 1 } });

  assert.equal(compareTruth(truth, untouched, { expect: "fresh" }).ok, true);
  const notTriaged = compareTruth(truth, untouched, { expect: "triaged" });
  assert.equal(notTriaged.ok, false);
  assert.deepEqual(notTriaged.samples[0].problems, ["TRUTH's bad frames are still kept: MO_H2O_S2"]);

  assert.equal(compareTruth(truth, triaged, { expect: "triaged" }).ok, true);
  const early = compareTruth(truth, triaged, { expect: "fresh" });
  assert.equal(early.ok, false);
  assert.deepEqual(early.samples[0].problems, ["frames dropped before triage: MO_H2O_S2"]);

  assert.equal(compareTruth(truth, untouched).ok, true);
  assert.equal(compareTruth(truth, triaged).ok, true);
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
    { sample: "MO_H2O", shown: "unreviewed", confirmed: "Pn3m" },
  ]);

  const shown = compareTruth(TRUTH, app({ h2o: h2oAssigned, h2oShown: { phase: "Pn3m" } }));
  assert.equal(shown.invariants.status_matches_assignment.ok, true);

  const wrong = compareTruth(TRUTH, app({ h2o: h2oAssigned, h2oShown: { phase: "Ia3d" } }));
  assert.equal(wrong.invariants.status_matches_assignment.ok, false);

  const unconfirmedShown = compareTruth(TRUTH, app({ pegShown: { phase: "Im3m" } }));
  assert.deepEqual(unconfirmedShown.invariants.status_matches_assignment.violations, [{ sample: "MO_PEG5", shown: "Im3m", confirmed: null }]);

  const coexistence = { state: "indexed", phases: ["Pn3m", "Im3m"] };
  const peg = compareTruth(TRUTH, app({ peg: [frame("MO_PEG5_S3", ["Im3m"], { assignment: coexistence })], pegShown: { phase: "Pn3m" } }));
  assert.equal(peg.invariants.status_matches_assignment.ok, true, "the corpus shows one phase of a coexistence call");

  const onDroppedFrames = [frame("MO_H2O_S1", ["Pn3m"], { assignment: assigned, ...dropped }), frame("MO_H2O_S2", ["Pn3m"])];
  const onDropped = compareTruth(TRUTH, app({ h2o: onDroppedFrames }));
  assert.deepEqual(onDropped.invariants.status_matches_assignment.violations, [
    { sample: "MO_H2O", shown: "unreviewed", confirmed: "Pn3m" },
  ], "a call made on a frame that was later dropped still belongs to the sample");
  const onDroppedShown = compareTruth(TRUTH, app({ h2o: onDroppedFrames, h2oShown: { phase: "Pn3m" } }));
  assert.equal(onDroppedShown.invariants.status_matches_assignment.ok, true);

  const acrossFrames = [frame("MO_H2O_S1", ["Pn3m"], { assignment: assigned }), frame("MO_H2O_S2", ["Im3m"], { assignment: { state: "indexed", phases: ["Im3m"] } })];
  const across = compareTruth(TRUTH, app({ h2o: acrossFrames, h2oShown: { phase: "Im3m" } }));
  assert.equal(across.invariants.status_matches_assignment.ok, true, "phases picked on different frames add up to one call");

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

function withViews(samples, views, corpusText = "Not indexed", workspaceText = "Unreviewed", corpusKept = {}) {
  return samples.map((s) => {
    const last = s.frames.at(-1).filename;
    return {
      ...s,
      views: views[s.id] ?? { workspace: last },
      corpusText,
      workspaceText,
      corpusKept: s.id in corpusKept ? corpusKept[s.id] : `${s.frames.filter((f) => f.status !== "rejected").length}/${s.frames.length}`,
    };
  });
}

test("with --ui, the corpus Kept column must show the app's kept / total", () => {
  const h2o = [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", ["Pn3m"], dropped)];
  const same = compareTruth(TRUTH, withViews(app({ h2o, h2oShown: { representative_id: 1 } }), { 1: { workspace: "MO_H2O_S1" } }));
  assert.equal(same.ok, true);
  assert.equal(same.samples[0].shown.corpus_kept, "1/2");

  const stale = compareTruth(TRUTH, withViews(app({ h2o, h2oShown: { representative_id: 1 } }), { 1: { workspace: "MO_H2O_S1" } }, "Not indexed", "Unreviewed", { 1: "2/2" }));
  assert.equal(stale.ok, false);
  assert.deepEqual(stale.samples[0].problems, ['the corpus Kept column reads "2/2" but the app has 1/2 frames kept']);

  const unread = compareTruth(TRUTH, withViews(app(), {}, "Not indexed", "Unreviewed", { 2: null }));
  assert.equal(unread.ok, false);
  assert.match(unread.samples[1].problems[0], /Kept column reads "nothing"/);
});

test("views_agree compares the corpus frame with the frame the workspace opens on", () => {
  const agree = compareTruth(TRUTH, withViews(app(), {}));
  assert.equal(agree.invariants.views_agree.checked, true);
  assert.equal(agree.invariants.views_agree.ok, true);

  const split = compareTruth(TRUTH, withViews(app(), { 1: { workspace: "MO_H2O_S1" } }));
  assert.equal(split.invariants.views_agree.ok, false);
  assert.deepEqual(split.invariants.views_agree.violations, [
    { sample: "MO_H2O", corpus_api: "MO_H2O_S2", workspace: "MO_H2O_S1" },
  ]);
  assert.deepEqual(split.samples[0].views, { workspace: "MO_H2O_S1" });
});

test("a view opening on a dropped frame breaks representative_kept", () => {
  const h2o = [frame("MO_H2O_S1", ["Pn3m"], dropped), frame("MO_H2O_S2", ["Pn3m"])];
  const r = compareTruth(TRUTH, withViews(app({ h2o }), { 1: { workspace: "MO_H2O_S1" } }));
  assert.deepEqual(r.invariants.representative_kept.violations, [
    { sample: "MO_H2O", view: "workspace", frame: "MO_H2O_S1", kept: ["MO_H2O_S2"] },
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

test("a workspace read that finds no frame is a problem, and the view invariants report unchecked", () => {
  const r = compareTruth(TRUTH, withViews(app(), { 1: { workspace: null }, 2: { workspace: "MO_PEG5_S3" } }));
  assert.equal(r.ok, false);
  assert.match(r.samples[0].problems[0], /could not read which frame the workspace opens on/);
  assert.deepEqual(r.samples[1].problems, []);
  for (const name of ["representative_kept", "views_agree", "workspace_status_matches"]) {
    assert.equal(r.invariants[name].checked, false);
    assert.equal(r.invariants[name].ok, null);
  }
  assert.deepEqual(r.invariants.views_agree.violations, []);
  assert.equal(r.invariants.status_matches_assignment.checked, true);
});

test("the workspace opening nothing is expected when every frame is dropped", () => {
  const h2o = [frame("MO_H2O_S1", ["Pn3m"], dropped), frame("MO_H2O_S2", ["Pn3m"], dropped)];
  const r = compareTruth(TRUTH, withViews(app({ h2o }), { 1: { workspace: null } }));
  assert.deepEqual(r.samples[0].problems, []);
  assert.equal(r.invariants.views_agree.checked, true);
});

test("with --ui, the workspace list must show the status GET /api/samples reports", () => {
  const same = compareTruth(TRUTH, withViews(app(), {}));
  assert.equal(same.invariants.workspace_status_matches.checked, true);
  assert.equal(same.invariants.workspace_status_matches.ok, true);
  assert.equal(same.samples[0].shown.workspace_text, "Unreviewed");

  const stale = compareTruth(TRUTH, withViews(app({ pegShown: { phase: "Im3m" } }), {}, "Im3m", "Unreviewed"));
  assert.equal(stale.ok, false);
  assert.ok(stale.summary.invariants_failed.includes("workspace_status_matches"));
  assert.deepEqual(stale.invariants.workspace_status_matches.violations.find((v) => v.sample === "MO_PEG5"), {
    sample: "MO_PEG5", shown: "Im3m", workspace: "Unreviewed",
  });
});

test("a workspace read that finds no current frame is a problem", () => {
  const r = compareTruth(TRUTH, withViews(app(), { 1: { workspace: null } }));
  assert.equal(r.ok, false);
  assert.match(r.samples[0].problems[0], /could not read which frame the workspace opens on/);
  assert.equal(r.invariants.workspace_status_matches.checked, false);
});

test("without --ui the workspace invariant is unchecked", () => {
  const r = compareTruth(TRUTH, app());
  assert.equal(r.invariants.workspace_status_matches.checked, false);
  assert.equal(r.invariants.workspace_status_matches.ok, null);
});
