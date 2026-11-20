import { getJson, request } from "./http.mjs";
import { findDemoExperiment, readTruth } from "./demo.mjs";
import { UsageError, intFlag } from "./args.mjs";
import { readViews, statusFromText } from "./views.mjs";

export function normalizeCall(phases) {
  const list = (Array.isArray(phases) ? phases : String(phases ?? "").split("+"))
    .map((p) => String(p).trim())
    .filter((p) => p && p !== "none");
  return list.length ? [...new Set(list)].sort().join("+") : "none";
}

export function isKept(frame) {
  return frame.status !== "rejected";
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

function statusMatches(shown, confirmed) {
  if (!confirmed.length) return shown.status === "unreviewed";
  if (shown.status === "unreviewed") return false;
  return confirmed.some((c) => (shown.status === "form_factor" ? c.call === "none" : c.call.split("+").includes(shown.call)));
}

const VIEW_NAMES = { loupe: "Loupe", focus: "Focus", workspace: "the workspace" };

export function unreadViews(app) {
  if (!app.views) return [];
  const unread = [];
  if (app.frames.length && !app.views.loupe) unread.push("loupe");
  if (app.frames.some(isKept) && !app.views.focus) unread.push("focus");
  if (app.frames.some(isKept) && !app.views.workspace) unread.push("workspace");
  return unread;
}

function sampleInvariants(app, frames) {
  const kept = frames.filter(isKept);
  const keptNames = kept.map((f) => f.filename);
  const shown = shownStatus(app.shown);
  const rep = frames.find((f) => f.id === (app.shown.representative_id ?? null)) ?? null;
  const confirmed = kept
    .map((f) => ({ frame: f.filename, call: assignedCall(f.assignment) }))
    .filter((c) => c.call !== null);
  const opened = { corpus_api: rep?.filename ?? null, ...(app.views ?? {}) };
  const onDropped = Object.entries(opened)
    .map(([view, name]) => ({ view, frame: frames.find((f) => f.filename === name) }))
    .filter(({ frame }) => frame && !isKept(frame) && kept.length)
    .map(({ view, frame }) => ({ view, frame: frame.filename, kept: keptNames }));
  const allRead = app.views && unreadViews(app).length === 0;
  const workspaceRead = app.views ? statusFromText(app.workspaceText) : null;
  return {
    workspace_status_matches: app.views && workspaceRead !== shown.status ? [{ shown: shownLabel(shown), workspace: app.workspaceText ?? null }] : null,
    status_matches_assignment: statusMatches(shown, confirmed) ? null : [{ shown: shownLabel(shown), confirmed }],
    representative_kept: onDropped.length ? onDropped : null,
    views_agree: allRead && kept.length && new Set(Object.values(opened)).size > 1 ? [opened] : null,
  };
}

export const INVARIANTS = {
  status_matches_assignment: "A sample with no confirmed call (no assignment on a kept frame) is shown as unreviewed; a sample with one is shown with that call.",
  representative_kept: "The frame the sample status describes (corpus_api) or Loupe, Focus or the workspace opens on (--ui) is never a dropped frame while the sample has a kept one.",
  views_agree: "The sample status (corpus_api), Loupe, Focus and the workspace describe the same frame of a sample (needs --ui).",
  workspace_status_matches: "The workspace sample list shows the same status as GET /api/samples (needs --ui).",
};

const NEEDS_VIEWS = new Set(["representative_kept", "views_agree", "workspace_status_matches"]);
const NEEDS_UI = new Set(["views_agree", "workspace_status_matches"]);

function rollupInvariants(perSample, ui) {
  const out = {};
  for (const [name, rule] of Object.entries(INVARIANTS)) {
    let checked = true;
    if (NEEDS_UI.has(name) && !ui.used) checked = false;
    if (NEEDS_VIEWS.has(name) && ui.unread) checked = false;
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

export function compareTruth(truthRows, appSamples, { strict = false, expect = null } = {}) {
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
    const autoCalls = frames.map((f) => (flared.has(f.filename)
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
      shown.corpus_text = one.corpusText;
      const read = statusFromText(one.corpusText);
      const expected = one.frames.length ? shown.status : "no_exposures";
      if (read !== expected) problems.push(`the corpus reads "${one.corpusText ?? "nothing"}" but GET /api/samples says ${expected}; check-truth is misreading the app`);
      shown.workspace_text = one.workspaceText ?? null;
      shown.corpus_kept = one.corpusKept ?? null;
      const keptCount = `${one.frames.filter(isKept).length}/${one.frames.length}`;
      if (one.corpusKept !== keptCount) problems.push(`the corpus Kept column reads "${one.corpusKept ?? "nothing"}" but the app has ${keptCount} frames kept`);
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
      invariants: one ? sampleInvariants(one, one.frames) : null,
    };
  });

  const failed = samples.filter((s) => !s.ok);
  const ui = { used: appSamples.some((s) => s.views), unread: appSamples.some((s) => unreadViews(s).length) };
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
  const [samples, exported] = await Promise.all([
    getJson(api, `/api/samples?experiment_id=${experimentId}`),
    getJson(api, `/api/experiments/${experimentId}/export`),
  ]);
  const autoByFrame = new Map();
  for (const s of exported) for (const e of s.exposures) autoByFrame.set(e.filename, e.indices.map((i) => i.phase));
  return Promise.all(
    samples.map(async (s) => {
      const exposures = await getJson(api, `/api/samples/${s.id}/exposures`);
      const frames = exposures
        .slice()
        .sort((a, b) => a.id - b.id)
        .map((e) => ({ id: e.id, filename: e.filename, status: e.status, selected: Boolean(e.selected), autoPhases: autoByFrame.get(e.filename) ?? [] }));
      for (const f of frames) f.assignment = await readAssignment(api, f.id);
      return {
        id: s.id,
        name: s.name,
        shown: { phase: s.phase ?? null, state: s.assignment_state ?? null, representative_id: s.representative_exposure_id ?? null },
        frameOrder: exposures.map((e) => e.filename),
        frames,
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
  let workspaceRead = null;
  if (flags.ui) {
    const opened = await views(cfg, experimentId, app.map((s) => ({ id: s.id, name: s.name, frameOrder: s.frameOrder })));
    const reused = [...opened.values()].map((v) => v.workspace_reused).filter((v) => v !== undefined);
    workspaceRead = reused.length ? (reused.every(Boolean) ? "in place" : "loaded") : null;
    for (const s of app) {
      const { corpus_text: corpusText = null, corpus_kept: corpusKept = null, workspace_text: workspaceText = null, workspace_reused: _reused, ...frames } = opened.get(s.id) ?? {};
      s.views = { loupe: frames.loupe ?? null, focus: frames.focus ?? null, workspace: frames.workspace ?? null };
      s.corpusText = corpusText;
      s.corpusKept = corpusKept;
      s.workspaceText = workspaceText;
    }
  }
  for (const s of app) delete s.frameOrder;
  const result = { experiment_id: experimentId, ...compareTruth(truth, app, { strict: Boolean(flags.strict), expect }) };
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
