import { getJson, request } from "./http.mjs";
import { findDemoExperiment, readTruth } from "./demo.mjs";
import { UsageError, intFlag } from "./args.mjs";
import { readViews, statusFromText } from "./views.mjs";
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
  if (shown.phase) return { status: "phase", call: normalizeCall([shown.phase]) };
  if (shown.state === "form_factor") return { status: "form_factor", call: "none" };
  return { status: "unreviewed", call: null };
}

function filenameOf(frames, id) {
  return frames.find((f) => f.id === id)?.filename ?? null;
}

function shownLabel(shown) {
  return shown.status === "phase" ? shown.call : shown.status;
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
  const expected = call === null ? "none" : normalizeCall(call.split("+").filter((p) => offered.has(p)));
  return !rail.form_factor && ticked === expected ? null : { rail: rail.form_factor ? "form_factor" : ticked, call: call ?? "unreviewed", expected };
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
  return shown.status === "phase" && call.split("+").includes(shown.call);
}

const VIEW_NAMES = { workspace: "the workspace" };

export function unreadViews(app) {
  if (!app.views) return [];
  const unread = [];
  if (app.frames.some(isKept) && !app.views.workspace) unread.push("workspace");
  return unread;
}

function sampleInvariants(app, frames, held = new Set()) {
  const kept = frames.filter(isKept);
  const keptNames = kept.map((f) => f.filename);
  const shown = shownStatus(app.shown);
  const rep = frames.find((f) => f.id === (app.shown.representative_id ?? null)) ?? null;
  const confirmed = sampleCall(frames);
  const opened = { corpus_api: rep?.filename ?? null, ...(app.views ?? {}) };
  const onDropped = Object.entries(opened)
    .map(([view, name]) => ({ view, frame: frames.find((f) => f.filename === name) }))
    .filter(({ frame }) => frame && !isKept(frame) && kept.length)
    .map(({ view, frame }) => ({ view, frame: frame.filename, kept: keptNames }));
  const allRead = app.views && unreadViews(app).length === 0;
  const workspaceRead = app.views ? statusFromText(app.workspaceText) : null;
  return {
    plot_matches_trace: app.views && app.workspacePlot !== undefined ? ((m) => (m ? [m] : null))(plotMismatch(app.views.workspace, app.workspacePlot, held)) : null,
    workspace_status_matches: app.views && workspaceRead !== shown.status ? [{ shown: shownLabel(shown), workspace: app.workspaceText ?? null }] : null,
    status_matches_assignment: statusMatches(shown, confirmed) ? null : [{ shown: shownLabel(shown), confirmed }],
    rail_matches_call: railViolation(app, kept, confirmed),
    representative_kept: onDropped.length ? onDropped : null,
    views_agree: allRead && kept.length && new Set(Object.values(opened)).size > 1 ? [opened] : null,
    bridged_labels_match: app.labelDrift?.length ? app.labelDrift : null,
    bridged_calls_match: app.callDrift?.length ? app.callDrift : null,
  };
}

export const INVARIANTS = {
  status_matches_assignment: "A sample's confirmed call is the phases assigned on any of its frames, kept or dropped, or the sample's form-factor state. A sample with none is shown as unreviewed; a sample with one is shown with that call.",
  rail_matches_call: "The workspace's Call list ticks exactly the sample's called phases that the opened frame offers as candidates, and marks form factor only for a form-factor call (needs --ui).",
  representative_kept: "The frame the sample status describes (corpus_api) or the workspace opens on (--ui) is never a dropped frame while the sample has a kept one.",
  views_agree: "The sample status (corpus_api) and the workspace describe the same frame of a sample (needs --ui).",
  workspace_status_matches: "The workspace sample list shows the same status as GET /api/samples (needs --ui).",
  bridged_labels_match: "Samples bridged as the same sample carry the same confirmed labels, apart from flight-path labels such as \"1.7 m\", which stay with each setup. Unchecked while the experiment has no bridge.",
  bridged_calls_match: "Samples bridged as the same sample show one call: the same call state and the same called phases. Unchecked while the experiment has no bridge.",
  plot_matches_trace: "The workspace plot of the opened frame says \"Waiting for this frame's trace.\" when its trace is held back in late/, and shows the trace otherwise (needs --ui).",
};

const NEEDS_VIEWS = new Set(["representative_kept", "views_agree", "workspace_status_matches", "rail_matches_call", "plot_matches_trace"]);
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

export const WAITING_PLOT = "Waiting for this frame's trace.";

export function plotMismatch(frame, plot, held) {
  if (!frame) return null;
  if (held.has(frame)) return plot === WAITING_PLOT ? null : { frame, plot, expected: WAITING_PLOT };
  return plot === null || plot === undefined ? null : { frame, plot, expected: "the trace" };
}

function noticeMismatch(held, notice) {
  if (notice === undefined) return null;
  const n = held.length;
  const says = /(\d+) frames? waiting for a trace/.exec(notice ?? "");
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
    const wrongAuto = autoCalls.filter((c) => (c.expected.startsWith("within ") ? !withinCall(c.call, expectedAuto) : c.call !== c.expected));
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
    const shown = one ? { ...shownStatus(one.shown), frame: filenameOf(one.frames, one.shown.representative_id ?? null) } : null;
    const unread = one ? unreadViews(one) : [];
    if (one && one.views) {
      shown.workspace_text = one.workspaceText ?? null;
      shown.workspace_kept = one.workspaceKept ?? null;
      const keptCount = one.frames.length ? `${one.frames.filter(isKept).length}/${one.frames.length}` : null;
      if (one.workspaceKept !== keptCount) problems.push(`the workspace's frame cards show ${one.workspaceKept ?? "nothing"} kept but the app has ${keptCount ?? "no frames"}`);
      for (const view of unread) problems.push(`could not read which frame ${VIEW_NAMES[view]} opens on (no single current "Frame N" thumbnail); check-truth is misreading the app`);
    }
    const assignedFrames = frames
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
  for (const s of exported) for (const e of s.exposures) {
    if (!Array.isArray(e.auto_phases)) throw new Error(`export frame ${e.filename} has no auto_phases array`);
    autoByFrame.set(e.filename, e.auto_phases);
  }
  return Promise.all(
    samples.map(async (s) => {
      const exposures = await getJson(api, `/api/samples/${s.id}/exposures`);
      const frames = exposures
        .slice()
        .sort((a, b) => a.id - b.id)
        .map((e) => ({ id: e.id, filename: e.filename, status: e.status, autoPhases: autoByFrame.get(e.filename) ?? [] }));
      for (const f of frames) f.assignment = await readAssignment(api, f.id);
      return {
        id: s.id,
        name: s.name,
        shown: { phase: s.phase ?? null, state: s.assignment_state ?? null, representative_id: s.representative_exposure_id ?? null },
        frameOrder: exposures.map((e) => e.filename),
        frames,
        bridged: bridged.has(s.id),
        ...(drift.has(s.id) ? { labelDrift: drift.get(s.id) } : {}),
        ...(callDrift.has(s.id) ? { callDrift: callDrift.get(s.id) } : {}),
      };
    }),
  );
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
      const { workspace_text: workspaceText = null, workspace_kept: workspaceKept = null, workspace_rail: workspaceRail = null, workspace_reused: _reused, workspace_plot: workspacePlot, workspace_notice: _notice, ...frames } = opened.get(s.id) ?? {};
      if (workspacePlot !== undefined) s.workspacePlot = workspacePlot;
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
