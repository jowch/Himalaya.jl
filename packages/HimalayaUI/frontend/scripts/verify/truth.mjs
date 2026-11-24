import { getJson, request } from "./http.mjs";
import { findDemoExperiment, readTruth } from "./demo.mjs";
import { UsageError, intFlag } from "./args.mjs";
import { callFromText, readViews, statusFromText } from "./views.mjs";
import { heldTraces } from "./late.mjs";

export function normalizeCall(phases) {
  const list = (Array.isArray(phases) ? phases : String(phases ?? "").split("+"))
    .map((p) => String(p).trim())
    .filter((p) => p && p !== "none");
  return list.length ? [...new Set(list)].sort().join("+") : "none";
}

export function isKept(frame) {
  return frame.status !== "rejected";
}

const DISTANCE_LABEL = /^\s*\d+(\.\d+)?\s*m\s*$/;

export function bridgeLabelDrift(samples, labels, bridges) {
  const shared = new Map(samples.map((s) => [s.id, new Set()]));
  for (const l of labels) {
    if (l.state === "confirmed" && !DISTANCE_LABEL.test(l.label)) shared.get(l.sample_id)?.add(l.label.trim().toLowerCase());
  }
  const groupOf = new Map(bridges.map((b) => [b.sample_id, b.group_id]));
  const nameOf = new Map(samples.map((s) => [s.id, s.name]));
  const key = (id) => [...(shared.get(id) ?? [])].sort();
  const out = new Map();
  for (const s of samples) {
    const g = groupOf.get(s.id);
    if (g === undefined) continue;
    const mine = key(s.id);
    const drift = bridges
      .filter((b) => b.group_id === g && b.sample_id !== s.id && nameOf.has(b.sample_id))
      .map((b) => ({ mate: nameOf.get(b.sample_id), mine, theirs: key(b.sample_id) }))
      .filter((d) => d.mine.join("\n") !== d.theirs.join("\n"));
    if (drift.length) out.set(s.id, drift);
  }
  return out;
}

export function bridgeCallDrift(samples, bridges) {
  const groupOf = new Map(bridges.map((b) => [b.sample_id, b.group_id]));
  const byId = new Map(samples.map((s) => [s.id, s]));
  const call = (s) => ({ state: s.assignment_state ?? "indexed", phases: [...(s.call_phases ?? [])].sort() });
  const same = (x, y) => x.state === y.state && x.phases.join("\n") === y.phases.join("\n");
  const out = new Map();
  for (const s of samples) {
    const g = groupOf.get(s.id);
    if (g === undefined) continue;
    const mine = call(s);
    const drift = bridges
      .filter((b) => b.group_id === g && b.sample_id !== s.id && byId.has(b.sample_id))
      .map((b) => ({ mate: byId.get(b.sample_id).name, mine, theirs: call(byId.get(b.sample_id)) }))
      .filter((d) => !same(d.mine, d.theirs));
    if (drift.length) out.set(s.id, drift);
  }
  return out;
}

function words(value) {
  return String(value ?? "").split(/\s+/).filter(Boolean);
}

function assignedCall(assignment) {
  if (!assignment) return null;
  if (assignment.state === "form_factor") return "none";
  if (assignment.state === "null") return null;
  return assignment.phases.length ? normalizeCall(assignment.phases) : null;
}

export function shownStatus(shown) {
  if (shown.phase) {
    const full = Array.isArray(shown.call_phases) && shown.call_phases.length > 0;
    return { status: "phase", call: normalizeCall(full ? shown.call_phases : [shown.phase]), ...(full ? { full } : {}) };
  }
  if (shown.state === "form_factor") return { status: "form_factor", call: "none" };
  return { status: "unreviewed", call: null };
}

function withoutFull({ full: _full, ...status }) {
  return status;
}

function filenameOf(frames, id) {
  return frames.find((f) => f.id === id)?.filename ?? null;
}

function shownLabel(shown) {
  return shown.status === "phase" ? shown.call : shown.status;
}

function callEvidenceOf(app) {
  return [...app.frames, ...(app.combined ? [app.combined] : []), ...(app.earlierAverages ?? [])];
}

export function withBridgedCallEvidence(apps, bridges) {
  const groupOf = new Map(bridges.map((b) => [b.sample_id, b.group_id]));
  return apps.map((app) => {
    const group = groupOf.get(app.id);
    const members = group === undefined ? [] : apps.filter((m) => groupOf.get(m.id) === group);
    return members.length > 1 ? { ...app, callEvidence: members.flatMap(callEvidenceOf) } : app;
  });
}

export function sampleCall(frames) {
  const states = new Set(frames.map((f) => f.assignment?.state).filter(Boolean));
  if (states.has("form_factor")) return "none";
  if (states.has("null")) return null;
  const phases = frames.flatMap((f) => f.assignment?.phases ?? []);
  return phases.length ? normalizeCall(phases) : null;
}

export function railMismatch(rail, call) {
  if (!rail) return null;
  const ticked = normalizeCall(rail.ticked);
  if (call === "none") return rail.form_factor && ticked === "none" ? null : { rail: rail.form_factor ? `form_factor+${ticked}` : ticked, call };
  const offered = new Set(rail.candidates);
  const flagged = new Set(rail.unsupported ?? []);
  const expected = call === null ? "none" : normalizeCall(call);
  const unmarked = rail.ticked.filter((p) => !offered.has(p) && !flagged.has(p));
  const marked = [...flagged].filter((p) => offered.has(p) || !rail.ticked.includes(p));
  if (!rail.form_factor && ticked === expected && !unmarked.length && !marked.length) return null;
  return {
    rail: rail.form_factor ? "form_factor" : ticked, call: call ?? "unreviewed", expected,
    ...(unmarked.length ? { unmarked } : {}), ...(marked.length ? { marked } : {}),
  };
}

function railViolation(app, kept, call) {
  if (!app.views || !kept.length) return null;
  if (!app.workspaceRail) return [{ rail: null, call: call ?? "unreviewed" }];
  const bad = railMismatch(app.workspaceRail, call);
  return bad ? [bad] : null;
}

function statusMatches(shown, call) {
  if (call === null) return shown.status === "unreviewed";
  if (shown.status === "form_factor") return call === "none";
  if (shown.status !== "phase") return false;
  return shown.full ? shown.call === call : call.split("+").includes(shown.call);
}

const VIEW_NAMES = { workspace: "the workspace" };

export function unreadViews(app) {
  if (!app.views) return [];
  const unread = [];
  if (app.frames.some(isKept) && !app.views.workspace) unread.push("workspace");
  return unread;
}

export const COMBINED = "combined";

export function isSelected(frame) {
  return isKept(frame) && frame.selected !== false;
}

export function selectedFrames(frames) {
  return frames.filter(isSelected);
}

export function selectionMismatch(app, frames, held = new Set()) {
  const selected = selectedFrames(frames);
  const names = (list) => list.map((f) => f.filename).sort().join(" ");
  const out = [];
  const analyzed = app.shown.analyzed_id ?? null;
  const dropped = frames.filter((f) => !isKept(f) && f.selected === true);
  if (dropped.length) out.push({ expected: "a dropped image is never selected", selected: names(dropped) });
  if (selected.length === 1 && analyzed !== selected[0].id) {
    out.push({ expected: `the sample status describes ${selected[0].filename}`, analyzed: frames.find((f) => f.id === analyzed)?.filename ?? analyzed });
  }
  const averaged = selected.length >= 2 && !selected.some((f) => held.has(f.filename));
  if (averaged) {
    if (!app.combined) out.push({ expected: `an average of ${names(selected)}`, combined: null });
    else {
      if (names(app.combined.sources) !== names(selected)) out.push({ expected: names(selected), combined: names(app.combined.sources) });
      if (analyzed !== app.combined.id) out.push({ expected: "the sample status describes the average", analyzed: frames.find((f) => f.id === analyzed)?.filename ?? analyzed });
    }
  }
  if (app.workspaceMembers !== undefined) {
    const drawn = [...app.workspaceMembers].sort((a, b) => a - b).join(" ");
    const want = averaged ? selected.map((f) => f.id).sort((a, b) => a - b).join(" ") : "";
    if (drawn !== want) out.push({ expected_plot: want ? `faint lines for ${names(selected)} under their average` : "one curve", drawn: app.workspaceMembers.map((id) => frames.find((f) => f.id === id)?.filename ?? id) });
  }
  return out.length ? out : null;
}

function sampleInvariants(app, frames, held = new Set()) {
  const kept = frames.filter(isKept);
  const keptNames = kept.map((f) => f.filename);
  const shown = shownStatus(app.shown);
  const repId = app.shown.analyzed_id ?? null;
  const rep = frames.find((f) => f.id === repId) ?? null;
  const confirmed = sampleCall(app.callEvidence ?? callEvidenceOf(app));
  const selected = selectedFrames(frames);
  const pending = selected.length >= 2 && selected.some((f) => held.has(f.filename));
  const noneSelected = kept.length > 0 && selected.length === 0;
  const opened = { corpus_api: app.combined && app.combined.id === repId ? COMBINED : rep?.filename ?? null, ...(app.views ?? {}) };
  const onDropped = Object.entries(opened)
    .map(([view, name]) => ({ view, frame: frames.find((f) => f.filename === name) }))
    .filter(({ frame }) => frame && !isKept(frame) && kept.length)
    .map(({ view, frame }) => ({ view, frame: frame.filename, kept: keptNames }));
  const allRead = app.views && unreadViews(app).length === 0;
  const workspaceRead = app.views ? statusFromText(app.workspaceText) : null;
  return {
    plot_matches_trace: app.views && app.workspacePlot !== undefined ? ((m) => (m ? [m] : null))(plotMismatch(app.views.workspace, app.workspacePlot, held, pending, noneSelected)) : null,
    workspace_status_matches: app.views && (workspaceRead !== shown.status || (shown.status === "phase" && callFromText(app.workspaceText) !== shown.call)) ? [{ shown: shownLabel(shown), workspace: app.workspaceText ?? null }] : null,
    status_matches_assignment: statusMatches(shown, confirmed) ? null : [{ shown: shownLabel(shown), confirmed }],
    rail_matches_call: railViolation(app, kept, confirmed),
    analyzed_kept: onDropped.length ? onDropped : null,
    views_agree: allRead && kept.length && !pending && !noneSelected && new Set(Object.values(opened)).size > 1 ? [opened] : null,
    bridged_labels_match: app.labelDrift?.length ? app.labelDrift : null,
    bridged_calls_match: app.callDrift?.length ? app.callDrift : null,
    plot_shows_selection: selectionMismatch(app, frames, held),
  };
}

export const INVARIANTS = {
  status_matches_assignment: "A sample's confirmed call is the phases assigned on any of its frames, kept or dropped, or the sample's form-factor state; a bridged sample's is read from the frames of every sample in its bridge. A sample with none is shown as unreviewed; a sample with one is shown with that call.",
  rail_matches_call: "The workspace's Assignment list ticks every one of the sample's called phases, marks as not found exactly the ticked phases the analyzed curve doesn't offer as candidates, and marks form factor only for a form-factor call (needs --ui).",
  analyzed_kept: "The image the sample status describes (corpus_api) or the workspace plots (--ui) is never a dropped image while the sample has a kept one.",
  views_agree: "The sample status (corpus_api) and the workspace describe the same curve of a sample: the average of its selected images, or its one selected image. Unchecked for a sample with kept images and none selected, which analyzes nothing (needs --ui).",
  workspace_status_matches: "The workspace sample list shows the same status as GET /api/samples (needs --ui).",
  bridged_labels_match: "Samples bridged as the same sample carry the same confirmed labels, apart from flight-path labels such as \"1.7 m\", which stay with each setup. Unchecked while the experiment has no bridge.",
  bridged_calls_match: "Samples bridged as the same sample show one call: the same call state and the same called phases. Unchecked while the experiment has no bridge.",
  plot_shows_selection: "A sample analyzes exactly its selected images, and a dropped image is never selected. With one selected, the sample status describes that image. With two or more, each with a trace, an average of exactly those images exists and the status describes it, and the workspace draws one faint line per selected image under it (with --ui).",
  plot_matches_trace: "The workspace plot of the opened frame says \"Waiting for this image's trace.\" when its trace is held back in late/, or \"Waiting for the selected images' traces.\" when it averages images one of which is held back, \"Select an image to analyze.\" when the sample has kept images and none selected, and shows the trace otherwise (needs --ui).",
};

const NEEDS_VIEWS = new Set(["analyzed_kept", "views_agree", "workspace_status_matches", "rail_matches_call", "plot_matches_trace"]);
const NEEDS_UI = new Set(["views_agree", "workspace_status_matches", "rail_matches_call", "plot_matches_trace"]);

function rollupInvariants(perSample, ui) {
  const out = {};
  for (const [name, rule] of Object.entries(INVARIANTS)) {
    let checked = true;
    if (NEEDS_UI.has(name) && !ui.used) checked = false;
    if (NEEDS_VIEWS.has(name) && ui.unread) checked = false;
    if ((name === "bridged_labels_match" || name === "bridged_calls_match") && !ui.bridged) checked = false;
    const violations = perSample
      .filter((s) => s.invariants && s.invariants[name])
      .flatMap((s) => s.invariants[name].map((v) => ({ sample: s.sample, ...v })));
    out[name] = { rule, checked, ok: checked ? violations.length === 0 : null, violations };
  }
  return out;
}

export const EXPECTATIONS = {
  fresh: "Before triage: every frame is kept, including TRUTH's bad frames.",
  triaged: "After triage: every frame in TRUTH's drop column is dropped.",
};

function withinCall(call, allowed) {
  if (call === "none") return true;
  const parts = new Set(allowed.split("+"));
  return call.split("+").every((p) => parts.has(p));
}

export const TRACES_RULE = "Exactly the frames whose trace is held back in late/ are reported as waiting for a trace, and no analysis failed. With --ui, the workspace notice counts exactly those frames.";

export const WAITING_PLOT = "Waiting for this image's trace.";

export const WAITING_AVERAGE = "Waiting for the selected images' traces.";

export const NONE_SELECTED_PLOT = "Select an image to analyze.";

export function plotMismatch(frame, plot, held, pending = false, noneSelected = false) {
  if (pending) return plot === WAITING_AVERAGE ? null : { frame, plot, expected: WAITING_AVERAGE };
  if (noneSelected) return plot === NONE_SELECTED_PLOT ? null : { frame, plot, expected: NONE_SELECTED_PLOT };
  if (!frame) return null;
  if (held.has(frame)) return plot === WAITING_PLOT ? null : { frame, plot, expected: WAITING_PLOT };
  return plot === null || plot === undefined ? null : { frame, plot, expected: "the trace" };
}

function noticeMismatch(held, notice) {
  if (notice === undefined) return null;
  const n = held.length;
  const says = /(\d+) images? waiting for a trace/.exec(notice ?? "");
  if (n === 0) return says ? notice : null;
  return says && Number(says[1]) === n ? null : notice ?? "no notice";
}

export function compareTraces(held, problems, notice) {
  const heldSet = new Set(held);
  const waiting = problems.filter((p) => p.problem === "no_trace").map((p) => p.filename).sort();
  const waitingSet = new Set(waiting);
  const failed = problems.filter((p) => p.problem !== "no_trace").map((p) => ({ frame: p.filename, reason: p.reason }));
  const stillWaiting = waiting.filter((f) => !heldSet.has(f));
  const notReported = [...heldSet].filter((f) => !waitingSet.has(f)).sort();
  const badNotice = noticeMismatch([...heldSet], notice);
  return {
    rule: TRACES_RULE,
    ok: stillWaiting.length === 0 && notReported.length === 0 && failed.length === 0 && badNotice === null,
    ...(notice !== undefined ? { notice: notice ?? null, notice_ok: badNotice === null } : {}),
    held: [...heldSet].sort(),
    waiting,
    still_waiting: stillWaiting,
    not_reported: notReported,
    failed,
  };
}

export function compareTruth(truthRows, appSamples, { strict = false, expect = null, held = new Set() } = {}) {
  const ownerOf = new Map();
  for (const row of truthRows) for (const base of words(row.bases)) ownerOf.set(base, row.sample);
  const appByTruth = new Map();
  const strays = [];
  for (const s of appSamples) {
    const owners = new Set(s.frames.map((f) => ownerOf.get(f.filename) ?? null));
    if (owners.size !== 1 || owners.has(null)) {
      strays.push({ sample_id: s.id, name: s.name, owners: [...owners] });
      continue;
    }
    const label = [...owners][0];
    if (!appByTruth.has(label)) appByTruth.set(label, []);
    appByTruth.get(label).push(s);
  }

  const samples = truthRows.map((row) => {
    const truth = normalizeCall(row.phases);
    const expectedAuto = normalizeCall(row.auto);
    const expectedFrames = words(row.bases);
    const bad = new Set(words(row.drop));
    const flared = new Set(words(row.flare));
    const empty = new Set([...bad].filter((b) => !flared.has(b)));
    const matches = appByTruth.get(row.sample) ?? [];
    const problems = [];
    if (matches.length === 0) problems.push("sample missing from the app");
    if (matches.length > 1) problems.push(`split into ${matches.length} app samples`);
    const frames = matches.flatMap((s) => s.frames);
    const seen = frames.map((f) => f.filename);
    const missing = expectedFrames.filter((b) => !seen.includes(b));
    const kept = frames.filter(isKept);
    if (matches.length && missing.length) problems.push(`missing frames: ${missing.join(" ")}`);
    const autoCalls = frames.map((f) => (held.has(f.filename)
      ? { frame: f.filename, call: normalizeCall(f.autoPhases), expected: "none", trace: "held" }
      : flared.has(f.filename)
      ? { frame: f.filename, call: normalizeCall(f.autoPhases), expected: `within ${expectedAuto}` }
      : { frame: f.filename, call: normalizeCall(f.autoPhases), expected: empty.has(f.filename) ? "none" : expectedAuto }));
    const solo = matches.length === 1 ? matches[0] : null;
    if (solo?.combined) {
      const sources = new Set(solo.combined.sources.map((f) => f.filename));
      const noisy = [...sources].some((f) => bad.has(f));
      const allowed = normalizeCall([...expectedAuto.split("+"), ...truth.split("+")]);
      autoCalls.push({ frame: COMBINED, call: normalizeCall(solo.combined.autoPhases), expected: noisy ? `within ${allowed}` : `${expectedAuto}, within ${allowed}`, allowed, needs: noisy ? "none" : expectedAuto });
    }
    const autoOk = (c) => {
      if (c.allowed) return withinCall(c.call, c.allowed) && (c.needs === "none" || withinCall(c.needs, c.call));
      return c.expected.startsWith("within ") ? withinCall(c.call, expectedAuto) : c.call === c.expected;
    };
    const wrongAuto = autoCalls.filter((c) => !autoOk(c));
    if (wrongAuto.length) {
      problems.push(`auto call ${wrongAuto.map((c) => `${c.frame}=${c.call} (expected ${c.expected})`).join(" ")}`);
    }
    const knownGap = expectedAuto !== truth;
    if (strict && knownGap) problems.push(`auto call ${expectedAuto} misses the generated truth ${truth}`);
    const badKept = kept.filter((f) => bad.has(f.filename)).map((f) => f.filename);
    if (expect === "triaged" && badKept.length) problems.push(`TRUTH's bad frames are still kept: ${badKept.join(" ")}`);
    const droppedEarly = frames.filter((f) => !isKept(f)).map((f) => f.filename);
    if (expect === "fresh" && droppedEarly.length) problems.push(`frames dropped before triage: ${droppedEarly.join(" ")}`);
    const one = matches.length === 1 ? matches[0] : null;
    const analyzedId = one?.shown.analyzed_id ?? null;
    const shown = one ? { ...withoutFull(shownStatus(one.shown)), frame: one.combined && one.combined.id === analyzedId ? COMBINED : filenameOf(one.frames, analyzedId), exposure_id: analyzedId } : null;
    if (one && one.shown.phase && !(Array.isArray(one.shown.call_phases) && one.shown.call_phases.length > 0)) {
      problems.push(`GET /api/samples gives ${one.name} the phase ${one.shown.phase} but no call_phases, so check-truth can't tell whether the list shows every phase of the call`);
    }
    const unread = one ? unreadViews(one) : [];
    if (one && one.views) {
      shown.workspace_text = one.workspaceText ?? null;
      shown.workspace_kept = one.workspaceKept ?? null;
      const keptCount = one.frames.length ? `${one.frames.filter(isKept).length}/${one.frames.length}` : null;
      if (one.workspaceKept !== keptCount) problems.push(`the workspace's frame cards show ${one.workspaceKept ?? "nothing"} kept but the app has ${keptCount ?? "no frames"}`);
      for (const view of unread) problems.push(`could not read which frame ${VIEW_NAMES[view]} opens on (no single current "Image N" thumbnail); check-truth is misreading the app`);
    }
    const averages = matches.flatMap((s) => [...(s.combined ? [s.combined] : []), ...(s.earlierAverages ?? [])]);
    const assignedFrames = [...frames, ...averages]
      .map((f) => ({ frame: f.filename, call: assignedCall(f.assignment) }))
      .filter((f) => f.call !== null);
    return {
      sample: row.sample,
      sample_id: one ? one.id : null,
      ok: problems.length === 0,
      problems,
      truth,
      expected_auto: expectedAuto,
      known_gap: knownGap ? row.note || `auto ${expectedAuto} vs truth ${truth}` : null,
      frames: {
        expected: expectedFrames.length,
        found: frames.length,
        kept: kept.length,
        missing,
        to_drop: [...bad],
        bad_still_kept: badKept,
        good_dropped: frames.filter((f) => !isKept(f) && !bad.has(f.filename)).map((f) => f.filename),
      },
      auto_calls: autoCalls,
      shown,
      shown_matches_truth: shown && shown.call !== null ? shown.call === truth : null,
      assigned_frames: assignedFrames,
      ...(one && one.views ? { views: one.views } : {}),
      invariants: one ? sampleInvariants(one, one.frames, held) : null,
    };
  });

  const failed = samples.filter((s) => !s.ok);
  const ui = { used: appSamples.some((s) => s.views), unread: appSamples.some((s) => unreadViews(s).length), bridged: appSamples.some((s) => s.bridged) };
  const invariants = rollupInvariants(samples, ui);
  const brokenInvariants = Object.entries(invariants).filter(([, v]) => v.ok === false).map(([k]) => k);
  return {
    ok: failed.length === 0 && strays.length === 0 && brokenInvariants.length === 0,
    gate: strict ? "truth" : "auto",
    expect: expect ? { [expect]: EXPECTATIONS[expect] } : null,
    summary: {
      samples: samples.length,
      passed: samples.length - failed.length,
      failed: failed.map((s) => s.sample),
      known_gaps: samples.filter((s) => s.known_gap).map((s) => s.sample),
      frames_expected: samples.reduce((n, s) => n + s.frames.expected, 0),
      frames_found: samples.reduce((n, s) => n + s.frames.found, 0),
      frames_kept: samples.reduce((n, s) => n + s.frames.kept, 0),
      bad_still_kept: samples.flatMap((s) => s.frames.bad_still_kept),
      shown_reviewed: samples.filter((s) => s.shown && s.shown.status !== "unreviewed").length,
      invariants_failed: brokenInvariants,
    },
    invariants,
    strays,
    samples: samples.map(({ invariants: _, ...s }) => s),
  };
}

async function readAssignment(api, exposureId) {
  const a = await getJson(api, `/api/exposures/${exposureId}/assignment`);
  if (!a.members.length) return { state: a.state, phases: [] };
  const indices = await getJson(api, `/api/exposures/${exposureId}/indices`);
  const byId = new Map(indices.map((i) => [i.id, i.phase]));
  return { state: a.state, phases: a.members.map((m) => byId.get(m)).filter(Boolean) };
}

export async function collectApp(cfg, experimentId) {
  const api = cfg.apiBase;
  const [samples, exported, labels, bridges] = await Promise.all([
    getJson(api, `/api/samples?experiment_id=${experimentId}`),
    getJson(api, `/api/experiments/${experimentId}/export`),
    getJson(api, `/api/experiments/${experimentId}/labels`),
    getJson(api, `/api/experiments/${experimentId}/bridges`),
  ]);
  const drift = bridgeLabelDrift(samples, labels, bridges);
  const callDrift = bridgeCallDrift(samples, bridges);
  const groupSize = new Map();
  for (const b of bridges) groupSize.set(b.group_id, (groupSize.get(b.group_id) ?? 0) + 1);
  const bridged = new Set(bridges.filter((b) => groupSize.get(b.group_id) > 1).map((b) => b.sample_id));
  const autoByFrame = new Map();
  const autoByCombination = new Map();
  for (const s of exported) for (const e of s.exposures) {
    if (!Array.isArray(e.auto_phases)) throw new Error(`export frame ${e.filename ?? e.kind} has no auto_phases array`);
    if (e.kind === COMBINED) autoByCombination.set([...(e.sources ?? [])].sort().join(" "), e.auto_phases);
    else autoByFrame.set(e.filename, e.auto_phases);
  }
  const apps = await Promise.all(
    samples.map(async (s) => {
      const all = await getJson(api, `/api/samples/${s.id}/exposures`);
      const exposures = all.filter((e) => e.kind !== COMBINED);
      const frames = exposures
        .slice()
        .sort((a, b) => a.id - b.id)
        .map((e) => ({ id: e.id, filename: e.filename, status: e.status, selected: e.status !== "rejected" && e.selected !== false, autoPhases: autoByFrame.get(e.filename) ?? [] }));
      for (const f of frames) f.assignment = await readAssignment(api, f.id);
      const row = all.find((e) => e.id === s.analyzed_exposure_id && e.kind === COMBINED);
      const sourceIds = new Set((row?.sources ?? []).map((x) => x.source_exposure_id));
      const sources = frames.filter((f) => sourceIds.has(f.id));
      const combined = row
        ? { id: row.id, filename: COMBINED, sources, autoPhases: autoByCombination.get(sources.map((f) => f.filename).sort().join(" ")) ?? [], assignment: await readAssignment(api, row.id) }
        : null;
      const earlierAverages = await Promise.all(all.filter((e) => e.kind === COMBINED && e.id !== row?.id)
        .map(async (e) => ({ id: e.id, filename: COMBINED, assignment: await readAssignment(api, e.id) })));
      return {
        id: s.id,
        name: s.name,
        shown: { phase: s.phase ?? null, call_phases: s.call_phases, state: s.assignment_state ?? null, analyzed_id: s.analyzed_exposure_id ?? null },
        frameOrder: exposures.map((e) => e.filename),
        frames,
        combined,
        earlierAverages,
        bridged: bridged.has(s.id),
        ...(drift.has(s.id) ? { labelDrift: drift.get(s.id) } : {}),
        ...(callDrift.has(s.id) ? { callDrift: callDrift.get(s.id) } : {}),
      };
    }),
  );
  return withBridgedCallEvidence(apps, bridges);
}

async function resolveExperiment(cfg, flags) {
  const experimentId = intFlag(flags, "experiment", null);
  if (experimentId !== null) return experimentId;
  const exp = await findDemoExperiment(cfg).catch(() => {
    throw new UsageError(`backend not answering at ${cfg.apiBase}`, "run `bin/verify seed`");
  });
  if (!exp) throw new UsageError("the demo experiment is not ingested", "run `bin/verify seed`");
  if (exp.ingest_status !== "complete") throw new UsageError(`demo experiment ${exp.id} is ${exp.ingest_status}`, "wait for the scan to finish; `bin/verify doctor` shows progress");
  return exp.id;
}

function loadTruth(cfg, flags) {
  const truth = readTruth(cfg, flags.truth);
  if (!truth) throw new UsageError(`no TRUTH.csv at ${flags.truth ?? cfg.demoDir}`, "run `make demo-data` or pass --truth <file>");
  return truth;
}

export async function checkTruth(cfg, flags, { views = readViews } = {}) {
  const truth = loadTruth(cfg, flags);
  const experimentId = await resolveExperiment(cfg, flags);
  const expect = flags.expect ?? null;
  if (expect !== null && !(expect in EXPECTATIONS)) throw new UsageError(`unknown --expect ${expect}`, `use one of ${Object.keys(EXPECTATIONS).join(", ")}`);
  const app = await collectApp(cfg, experimentId);
  const held = heldTraces(cfg);
  const heldSet = new Set(held);
  let workspaceRead = null;
  let opened = null;
  if (flags.ui) {
    const heldIds = app.flatMap((s) => s.frames).filter((f) => heldSet.has(f.filename)).map((f) => f.id);
    opened = await views(cfg, experimentId, app.map((s) => ({ id: s.id, name: s.name, frameOrder: s.frameOrder, heldIds })), { reload: Boolean(flags.reload) });
    const reused = [...opened.values()].map((v) => v.workspace_reused).filter((v) => v !== undefined);
    workspaceRead = reused.length ? (reused.every(Boolean) ? "in place" : "loaded") : null;
    for (const s of app) {
      const { workspace_text: workspaceText = null, workspace_kept: workspaceKept = null, workspace_rail: workspaceRail = null, workspace_reused: _reused, workspace_plot: workspacePlot, workspace_notice: _notice, workspace_members: workspaceMembers, ...frames } = opened.get(s.id) ?? {};
      if (workspacePlot !== undefined) s.workspacePlot = workspacePlot;
      if (workspaceMembers !== undefined) s.workspaceMembers = workspaceMembers;
      s.views = { workspace: frames.workspace ?? null };
      s.workspaceKept = workspaceKept;
      s.workspaceText = workspaceText;
      s.workspaceRail = workspaceRail;
    }
  }
  for (const s of app) delete s.frameOrder;
  const problems = await getJson(cfg.apiBase, `/api/experiments/${experimentId}/problems`);
  const notices = flags.ui ? [...(opened?.values() ?? [])].map((v) => v.workspace_notice).filter((v) => v !== undefined) : [];
  const traces = compareTraces(held, problems.exposures ?? [], notices.length ? notices[notices.length - 1] : undefined);
  const compared = compareTruth(truth, app, { strict: Boolean(flags.strict), expect, held: heldSet });
  const result = { experiment_id: experimentId, ...compared, ok: compared.ok && traces.ok, traces };
  if (flags.ui) result.workspace_read = workspaceRead;
  return result;
}

export async function triage(cfg, flags) {
  const truth = loadTruth(cfg, flags);
  const experimentId = await resolveExperiment(cfg, flags);
  const app = await collectApp(cfg, experimentId);
  const bad = new Set(truth.flatMap((r) => words(r.drop)));
  const status = flags.undo ? null : "rejected";
  const changed = [];
  for (const s of app) {
    for (const f of s.frames) {
      if (!bad.has(f.filename) || (f.status === "rejected") === (status === "rejected")) continue;
      if (!flags["dry-run"]) {
        const res = await request(cfg.apiBase, "PATCH", `/api/exposures/${f.id}/status`, { data: { status } });
        if (res.status >= 300) throw new Error(`PATCH /api/exposures/${f.id}/status returned ${res.status}`);
      }
      changed.push({ sample_id: s.id, frame: f.filename, exposure_id: f.id, status: status ?? "kept" });
    }
  }
  const unknown = [...bad].filter((b) => !app.some((s) => s.frames.some((f) => f.filename === b)));
  return { ok: unknown.length === 0, experiment_id: experimentId, dry_run: Boolean(flags["dry-run"]), changed, ...(unknown.length ? { error: `frames in TRUTH's drop column are not in the app: ${unknown.join(" ")}` } : {}) };
}
