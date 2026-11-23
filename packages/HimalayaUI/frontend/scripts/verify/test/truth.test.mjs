import { test } from "node:test";
import assert from "node:assert/strict";
import { checkTruth, bridgeCallDrift, bridgeLabelDrift, compareTruth, selectedFrames, normalizeCall, railMismatch, shownStatus } from "../truth.mjs";

const TRUTH = [
  { sample: "MO_H2O", bases: "MO_H2O_S1 MO_H2O_S2", phases: "Pn3m", auto: "Pn3m", drop: "", note: "" },
  { sample: "MO_PEG5", bases: "MO_PEG5_S3", phases: "Pn3m+Im3m", auto: "Im3m", drop: "", note: "Pn3m √6 overlaps Im3m √10" },
];

const NOT_INDEXED = { phase: null, state: "indexed" };

function frame(filename, autoPhases, extra = {}) {
  return { id: Number(filename.at(-1)), filename, status: null, autoPhases, assignment: { state: "indexed", phases: [] }, ...extra };
}

function combinedOf(id, frames) {
  const sources = selectedFrames(frames);
  if (sources.length < 2) return null;
  const autoPhases = sources.find((f) => f.autoPhases.length)?.autoPhases ?? [];
  return { id: 100 + id, filename: "combined", sources, autoPhases, assignment: { state: "indexed", phases: [] } };
}

function sample(id, name, frames, shown = {}, combined = combinedOf(id, frames)) {
  const one = selectedFrames(frames);
  const rep = combined ? combined.id : one.length === 1 ? one[0].id : frames.at(-1).id;
  return { id, name, frames, combined, shown: { ...NOT_INDEXED, analyzed_id: rep, ...shown } };
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
  const triaged = app({ h2o: [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", [], dropped)], h2oShown: { analyzed_id: 1 } });

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
  assert.deepEqual(r.samples[1].shown, { status: "phase", call: "Im3m", frame: "MO_PEG5_S3", exposure_id: 3 });
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
  assert.equal(peg.ok, false, "a payload without call_phases can't show whether the list names every phase");
  assert.match(peg.samples[1].problems.join("\n"), /GET \/api\/samples gives MO_PEG5_S3 \(S01P02\) the phase Pn3m but no call_phases/);
  const pegAll = compareTruth(TRUTH, app({ peg: [frame("MO_PEG5_S3", ["Im3m"], { assignment: coexistence })], pegShown: { phase: "Pn3m", call_phases: ["Pn3m", "Im3m"] } }));
  assert.equal(pegAll.invariants.status_matches_assignment.ok, true);
  const pegShort = compareTruth(TRUTH, app({ peg: [frame("MO_PEG5_S3", ["Im3m"], { assignment: coexistence })], pegShown: { phase: "Pn3m", call_phases: ["Pn3m"] } }));
  assert.deepEqual(pegShort.invariants.status_matches_assignment.violations, [{ sample: "MO_PEG5", shown: "Pn3m", confirmed: "Im3m+Pn3m" }]);

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

test("analyzed_kept: the corpus status describes a dropped frame while a kept one exists", () => {
  const h2o = [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", ["Pn3m"], dropped)];
  const r = compareTruth(TRUTH, app({ h2o, h2oShown: { analyzed_id: 2 } }));
  assert.deepEqual(r.invariants.analyzed_kept.violations, [
    { sample: "MO_H2O", view: "corpus_api", frame: "MO_H2O_S2", kept: ["MO_H2O_S1"] },
  ]);

  const allDropped = compareTruth(TRUTH, app({ h2o: h2o.map((f) => ({ ...f, ...dropped })) }));
  assert.equal(allDropped.invariants.analyzed_kept.ok, true);

  const onKept = compareTruth(TRUTH, app({ h2o, h2oShown: { analyzed_id: 1 } }));
  assert.equal(onKept.invariants.analyzed_kept.ok, true);
});

function withViews(samples, views, workspaceText = "Unreviewed", workspaceKept = {}, rails = {}) {
  return samples.map((s) => {
    const last = s.combined ? "combined" : s.frames.at(-1).filename;
    return {
      ...s,
      views: views[s.id] ?? { workspace: last },
      workspaceText,
      workspaceRail: s.id in rails ? rails[s.id] : { ticked: [], candidates: [], form_factor: false },
      workspaceKept: s.id in workspaceKept ? workspaceKept[s.id] : `${s.frames.filter((f) => f.status !== "rejected").length}/${s.frames.length}`,
    };
  });
}

test("with --ui, the workspace's frame cards must show the app's kept / total", () => {
  const h2o = [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", ["Pn3m"], dropped)];
  const same = compareTruth(TRUTH, withViews(app({ h2o, h2oShown: { analyzed_id: 1 } }), { 1: { workspace: "MO_H2O_S1" } }));
  assert.equal(same.ok, true);
  assert.equal(same.samples[0].shown.workspace_kept, "1/2");

  const stale = compareTruth(TRUTH, withViews(app({ h2o, h2oShown: { analyzed_id: 1 } }), { 1: { workspace: "MO_H2O_S1" } }, "Unreviewed", { 1: "2/2" }));
  assert.equal(stale.ok, false);
  assert.deepEqual(stale.samples[0].problems, ["the workspace's frame cards show 2/2 kept but the app has 1/2"]);

  const unread = compareTruth(TRUTH, withViews(app(), {}, "Unreviewed", { 2: null }));
  assert.equal(unread.ok, false);
  assert.match(unread.samples[1].problems[0], /frame cards show nothing kept/);
});

test("views_agree compares the API status frame with the frame the workspace opens on", () => {
  const agree = compareTruth(TRUTH, withViews(app(), {}));
  assert.equal(agree.invariants.views_agree.checked, true);
  assert.equal(agree.invariants.views_agree.ok, true);

  const split = compareTruth(TRUTH, withViews(app(), { 1: { workspace: "MO_H2O_S1" } }));
  assert.equal(split.invariants.views_agree.ok, false);
  assert.deepEqual(split.invariants.views_agree.violations, [
    { sample: "MO_H2O", corpus_api: "combined", workspace: "MO_H2O_S1" },
  ]);
  assert.deepEqual(split.samples[0].views, { workspace: "MO_H2O_S1" });
});

test("a view opening on a dropped frame breaks analyzed_kept", () => {
  const h2o = [frame("MO_H2O_S1", ["Pn3m"], dropped), frame("MO_H2O_S2", ["Pn3m"])];
  const r = compareTruth(TRUTH, withViews(app({ h2o }), { 1: { workspace: "MO_H2O_S1" } }));
  assert.deepEqual(r.invariants.analyzed_kept.violations, [
    { sample: "MO_H2O", view: "workspace", frame: "MO_H2O_S1", kept: ["MO_H2O_S2"] },
  ]);
});

test("a workspace read that finds no frame is a problem, and the view invariants report unchecked", () => {
  const r = compareTruth(TRUTH, withViews(app(), { 1: { workspace: null }, 2: { workspace: "MO_PEG5_S3" } }));
  assert.equal(r.ok, false);
  assert.match(r.samples[0].problems[0], /could not read which frame the workspace opens on/);
  assert.deepEqual(r.samples[1].problems, []);
  for (const name of ["analyzed_kept", "views_agree", "workspace_status_matches"]) {
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

  const stale = compareTruth(TRUTH, withViews(app({ pegShown: { phase: "Im3m" } }), {}, "Unreviewed"));
  assert.equal(stale.ok, false);
  assert.ok(stale.summary.invariants_failed.includes("workspace_status_matches"));
  assert.deepEqual(stale.invariants.workspace_status_matches.violations.find((v) => v.sample === "MO_PEG5"), {
    sample: "MO_PEG5", shown: "Im3m", workspace: "Unreviewed",
  });
});

test("with --ui, the workspace list must show every phase of the call, not just the dominant one", () => {
  const coexistence = { state: "indexed", phases: ["Pn3m", "Im3m"] };
  const peg = [frame("MO_PEG5_S3", ["Im3m"], { assignment: coexistence })];
  const pegShown = { phase: "Im3m", call_phases: ["Im3m", "Pn3m"] };
  const one = compareTruth(TRUTH, withViews(app({ peg, pegShown }), {}, "Im3m").filter((s) => s.id === 2));
  assert.deepEqual(one.invariants.workspace_status_matches.violations, [{ sample: "MO_PEG5", shown: "Im3m+Pn3m", workspace: "Im3m" }]);
  const both = compareTruth(TRUTH, withViews(app({ peg, pegShown }), {}, "Im3m + Pn3m").filter((s) => s.id === 2));
  assert.equal(both.invariants.workspace_status_matches.ok, true);
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

test("railMismatch wants every called phase ticked, and form factor only for a form-factor call", () => {
  const rail = (ticked, candidates = ["Pn3m", "Im3m"], ff = false, unsupported = []) => ({ ticked, candidates, unsupported, form_factor: ff });
  assert.equal(railMismatch(rail(["Pn3m"]), "Pn3m"), null);
  assert.deepEqual(railMismatch(rail(["Pn3m"]), "Ia3d+Pn3m"), { rail: "Pn3m", call: "Ia3d+Pn3m", expected: "Ia3d+Pn3m" });
  assert.deepEqual(railMismatch(rail([]), "Ia3d"), { rail: "none", call: "Ia3d", expected: "Ia3d" });
  assert.deepEqual(railMismatch(rail([]), "Pn3m"), { rail: "none", call: "Pn3m", expected: "Pn3m" });
  assert.deepEqual(railMismatch(rail(["Im3m"]), null), { rail: "Im3m", call: "unreviewed", expected: "none" });
  assert.equal(railMismatch(rail([], undefined, true), "none"), null);
  assert.deepEqual(railMismatch(rail([], undefined, false), "none"), { rail: "none", call: "none" });
  assert.deepEqual(railMismatch(rail([], undefined, true), null), { rail: "form_factor", call: "unreviewed", expected: "none" });
});

test("railMismatch wants a called phase the curve doesn't offer ticked and marked, and only that one marked", () => {
  const rail = (ticked, unsupported, candidates = ["Pn3m", "Im3m"]) => ({ ticked, candidates, unsupported, form_factor: false });
  assert.equal(railMismatch(rail(["Ia3d", "Pn3m"], ["Ia3d"]), "Ia3d+Pn3m"), null);
  assert.deepEqual(railMismatch(rail(["Ia3d", "Pn3m"], []), "Ia3d+Pn3m"),
    { rail: "Ia3d+Pn3m", call: "Ia3d+Pn3m", expected: "Ia3d+Pn3m", unmarked: ["Ia3d"] });
  assert.deepEqual(railMismatch(rail(["Pn3m"], ["Pn3m"]), "Pn3m"),
    { rail: "Pn3m", call: "Pn3m", expected: "Pn3m", marked: ["Pn3m"] });
  assert.deepEqual(railMismatch(rail([], ["Ia3d"]), null),
    { rail: "none", call: "unreviewed", expected: "none", marked: ["Ia3d"] });
});

test("with --ui, a Call list that misses a phase called on a dropped frame breaks rail_matches_call", () => {
  const h2o = [frame("MO_H2O_S1", ["Pn3m"], { ...dropped, assignment: { state: "indexed", phases: ["Pn3m"] } }), frame("MO_H2O_S2", ["Pn3m"])];
  const shown = { h2o, h2oShown: { phase: "Pn3m", state: "indexed", analyzed_id: 2 } };
  const ok = compareTruth(TRUTH, withViews(app(shown), {}, "Pn3m", {}, { 1: { ticked: ["Pn3m"], candidates: ["Pn3m"], form_factor: false } }));
  assert.equal(ok.invariants.rail_matches_call.checked, true);
  assert.equal(ok.invariants.rail_matches_call.ok, true);

  const missed = compareTruth(TRUTH, withViews(app(shown), {}, "Pn3m", {}, { 1: { ticked: [], candidates: ["Pn3m"], form_factor: false } }));
  assert.equal(missed.invariants.rail_matches_call.ok, false);
  assert.deepEqual(missed.invariants.rail_matches_call.violations, [{ sample: "MO_H2O", rail: "none", call: "Pn3m", expected: "Pn3m" }]);

  const unread = compareTruth(TRUTH, withViews(app(shown), {}, "Pn3m", {}, { 1: null }));
  assert.equal(unread.invariants.rail_matches_call.ok, false);
});

test("bridgeLabelDrift compares confirmed labels across a bridge, ignoring case and flight paths", () => {
  const samples = [{ id: 1, name: "A" }, { id: 2, name: "B" }, { id: 3, name: "C" }];
  const bridges = [{ sample_id: 1, group_id: 7 }, { sample_id: 2, group_id: 7 }];
  const l = (sample_id, label, state = "confirmed") => ({ sample_id, label, state });
  const same = [l(1, "Buffer"), l(2, "buffer"), l(1, "0.5 m"), l(2, "1.7 m"), l(1, "maybe", "suggested"), l(3, "37C")];
  assert.equal(bridgeLabelDrift(samples, same, bridges).size, 0);
  const apart = bridgeLabelDrift(samples, [l(1, "Buffer"), l(3, "Buffer")], bridges);
  assert.deepEqual([...apart.keys()], [1, 2]);
  assert.deepEqual(apart.get(2), [{ mate: "A", mine: [], theirs: ["buffer"] }]);
});

test("bridged_labels_match is unchecked without a bridge, and fails on a sample whose labels drifted from its bridge", () => {
  const none = compareTruth(TRUTH, app());
  assert.equal(none.invariants.bridged_labels_match.checked, false);
  assert.equal(none.invariants.bridged_labels_match.ok, null);
  const bridged = app();
  bridged[0].bridged = true;
  const ok = compareTruth(TRUTH, bridged);
  assert.equal(ok.invariants.bridged_labels_match.checked, true);
  assert.equal(ok.invariants.bridged_labels_match.ok, true);
  const drifted = app();
  drifted[0].bridged = true;
  drifted[0].labelDrift = [{ mate: "MO_PEG5_S3 (S01P02)", mine: ["buffer"], theirs: [] }];
  const r = compareTruth(TRUTH, drifted);
  assert.deepEqual(r.summary.invariants_failed, ["bridged_labels_match"]);
  assert.deepEqual(r.invariants.bridged_labels_match.violations, [
    { sample: "MO_H2O", mate: "MO_PEG5_S3 (S01P02)", mine: ["buffer"], theirs: [] },
  ]);
});

test("bridgeCallDrift compares the call state and called phases across a bridge", () => {
  const sample = (id, name, assignment_state, call_phases) => ({ id, name, assignment_state, call_phases });
  const bridges = [{ sample_id: 1, group_id: 7 }, { sample_id: 2, group_id: 7 }];
  const same = [sample(1, "A", "indexed", ["Pn3m", "Im3m"]), sample(2, "B", "indexed", ["Im3m", "Pn3m"]), sample(3, "C", "form_factor", [])];
  assert.equal(bridgeCallDrift(same, bridges).size, 0);
  const apart = bridgeCallDrift([sample(1, "A", "indexed", ["Pn3m"]), sample(2, "B", "form_factor", [])], bridges);
  assert.deepEqual([...apart.keys()], [1, 2]);
  assert.deepEqual(apart.get(1), [{ mate: "B", mine: { state: "indexed", phases: ["Pn3m"] }, theirs: { state: "form_factor", phases: [] } }]);
});

test("bridged_calls_match is unchecked without a bridge and fails on a sample whose call drifted from its bridge", () => {
  assert.equal(compareTruth(TRUTH, app()).invariants.bridged_calls_match.checked, false);
  const drifted = app();
  drifted[0].bridged = true;
  drifted[0].callDrift = [{ mate: "B", mine: { state: "indexed", phases: ["Pn3m"] }, theirs: { state: "form_factor", phases: [] } }];
  const r = compareTruth(TRUTH, drifted);
  assert.deepEqual(r.summary.invariants_failed, ["bridged_calls_match"]);
});

test("plot_shows_selection: two selected images need their average, and the status describes it", () => {
  const ok = compareTruth(TRUTH, app());
  assert.equal(ok.invariants.plot_shows_selection.ok, true);
  assert.equal(ok.invariants.plot_shows_selection.checked, true);

  const h2o = [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", ["Pn3m"])];
  const missing = compareTruth(TRUTH, [sample(1, "MO_H2O_S1 (S01P01)", h2o, {}, null), app()[1]]);
  assert.deepEqual(missing.invariants.plot_shows_selection.violations, [
    { sample: "MO_H2O", expected: "an average of MO_H2O_S1 MO_H2O_S2", combined: null },
  ]);

  const notShown = compareTruth(TRUTH, app({ h2oShown: { analyzed_id: 2 } }));
  assert.deepEqual(notShown.invariants.plot_shows_selection.violations, [
    { sample: "MO_H2O", expected: "the sample status describes the average", analyzed: "MO_H2O_S2" },
  ]);

  const stale = [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", ["Pn3m"], dropped), frame("MO_H2O_S3", ["Pn3m"])];
  const staleApp = app({ h2o: stale });
  staleApp[0].combined = { ...staleApp[0].combined, sources: stale };
  const r = compareTruth([{ ...TRUTH[0], bases: "MO_H2O_S1 MO_H2O_S2 MO_H2O_S3" }, TRUTH[1]], staleApp);
  assert.deepEqual(r.invariants.plot_shows_selection.violations, [
    { sample: "MO_H2O", expected: "MO_H2O_S1 MO_H2O_S3", combined: "MO_H2O_S1 MO_H2O_S2 MO_H2O_S3" },
  ]);
});

test("plot_shows_selection: one selected image is analyzed as itself, and a dropped image is never selected", () => {
  const one = [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", ["Pn3m"], { selected: false })];
  const r = compareTruth(TRUTH, app({ h2o: one }));
  assert.equal(r.samples[0].shown.frame, "MO_H2O_S1");
  assert.equal(r.invariants.plot_shows_selection.ok, true);
  assert.equal(app({ h2o: one })[0].combined, null);

  const wrong = compareTruth(TRUTH, app({ h2o: one, h2oShown: { analyzed_id: 2 } }));
  assert.deepEqual(wrong.invariants.plot_shows_selection.violations, [
    { sample: "MO_H2O", expected: "the sample status describes MO_H2O_S1", analyzed: "MO_H2O_S2" },
  ]);

  const bad = [frame("MO_H2O_S1", ["Pn3m"]), frame("MO_H2O_S2", ["Pn3m"], { ...dropped, selected: true })];
  assert.deepEqual(compareTruth(TRUTH, app({ h2o: bad })).invariants.plot_shows_selection.violations, [
    { sample: "MO_H2O", expected: "a dropped image is never selected", selected: "MO_H2O_S2" },
  ]);
});

test("plot_shows_selection: with --ui the workspace draws one faint line per selected image", () => {
  const drawn = withViews(app(), {});
  drawn[0].workspaceMembers = [1, 2];
  drawn[1].workspaceMembers = [];
  assert.equal(compareTruth(TRUTH, drawn).invariants.plot_shows_selection.ok, true);

  const single = withViews(app(), {});
  single[0].workspaceMembers = [];
  single[1].workspaceMembers = [];
  assert.deepEqual(compareTruth(TRUTH, single).invariants.plot_shows_selection.violations, [
    { sample: "MO_H2O", expected_plot: "faint lines for MO_H2O_S1 MO_H2O_S2 under their average", drawn: [] },
  ]);
});

test("plot_shows_selection and views_agree wait while a selected image's trace is held back", () => {
  const pending = withViews(app(), {});
  pending[0].combined = null;
  pending[0].shown = { ...pending[0].shown, analyzed_id: null };
  pending[0].views = { workspace: "MO_H2O_S1" };
  pending[0].workspaceMembers = [];
  pending[0].workspacePlot = "Waiting for the selected images' traces.";
  pending[1].workspaceMembers = [];
  const r = compareTruth(TRUTH, pending, { held: new Set(["MO_H2O_S2"]) });
  assert.equal(r.invariants.plot_shows_selection.ok, true);
  assert.equal(r.invariants.views_agree.ok, true);
  assert.deepEqual(r.invariants.plot_matches_trace.violations, []);
});

test("a call confirmed on an earlier average still counts as the sample's call", () => {
  const earlier = app({ h2oShown: { phase: "Pn3m", state: "indexed" } });
  earlier[0].earlierAverages = [{ id: 99, filename: "combined", assignment: { state: "indexed", phases: ["Pn3m"] } }];
  assert.equal(compareTruth(TRUTH, earlier).invariants.status_matches_assignment.ok, true);
  const without = app({ h2oShown: { phase: "Pn3m", state: "indexed" } });
  assert.equal(compareTruth(TRUTH, without).invariants.status_matches_assignment.ok, false);
});

test("the average's auto call may add a TRUTH phase but must keep the auto call", () => {
  const peg = [frame("MO_PEG5_S3", ["Im3m"]), frame("MO_PEG5_S4", ["Im3m"])];
  const truth = [TRUTH[0], { ...TRUTH[1], bases: "MO_PEG5_S3 MO_PEG5_S4" }];
  const better = app({ peg });
  better[1].combined = { ...better[1].combined, autoPhases: ["Im3m", "Pn3m"] };
  assert.equal(compareTruth(truth, better).ok, true);

  const lost = app({ peg });
  lost[1].combined = { ...lost[1].combined, autoPhases: [] };
  const r = compareTruth(truth, lost);
  assert.equal(r.ok, false);
  assert.match(r.samples[1].problems[0], /combined=none/);

  const extra = app({ peg });
  extra[1].combined = { ...extra[1].combined, autoPhases: ["Im3m", "Ia3d"] };
  assert.equal(compareTruth(truth, extra).ok, false);
});

test("--gate invariants reports TRUTH problems without failing on them, and still fails on a broken invariant", () => {
  const offTruth = app({ h2o: [frame("MO_H2O_S1", []), frame("MO_H2O_S2", [])] });
  const auto = compareTruth(TRUTH, offTruth);
  assert.equal(auto.ok, false);
  const gated = compareTruth(TRUTH, offTruth, { gate: "invariants" });
  assert.equal(gated.gate, "invariants");
  assert.deepEqual(gated.summary.failed, auto.summary.failed);
  assert.equal(gated.ok, true);
  const broken = compareTruth(TRUTH, app({ pegShown: { phase: "Im3m" } }), { gate: "invariants" });
  assert.equal(broken.ok, false);
  assert.deepEqual(broken.summary.invariants_failed, ["status_matches_assignment"]);
});

test("check-truth refuses an unknown gate and --gate invariants with --strict", async () => {
  await assert.rejects(checkTruth({}, { gate: "truth" }), /unknown --gate truth/);
  await assert.rejects(checkTruth({}, { gate: "invariants", strict: true }), /disagree/);
  await assert.rejects(checkTruth({}, { gate: "invariants", expect: "triaged" }), /--expect disagree/);
});
