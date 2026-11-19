import { getJson } from "./http.mjs";
import { findDemoExperiment, readTruth } from "./demo.mjs";
import { UsageError } from "./args.mjs";

export function normalizeCall(phases) {
  const list = (Array.isArray(phases) ? phases : String(phases ?? "").split("+"))
    .map((p) => String(p).trim())
    .filter((p) => p && p !== "none");
  return list.length ? [...new Set(list)].sort().join("+") : "none";
}

export function isKept(frame) {
  return frame.status !== "rejected";
}

export function representativeFrame(frames) {
  const kept = frames.filter(isKept);
  const selected = kept.filter((f) => f.selected);
  if (selected.length) return selected[selected.length - 1];
  return kept.length ? kept[kept.length - 1] : null;
}

function assignedCall(assignment) {
  if (!assignment) return null;
  if (assignment.state === "form_factor" || assignment.state === "null") return "none";
  return assignment.phases.length ? normalizeCall(assignment.phases) : null;
}

export function compareTruth(truthRows, appSamples, { strict = false } = {}) {
  const ownerOf = new Map();
  for (const row of truthRows) for (const base of row.bases.split(/\s+/).filter(Boolean)) ownerOf.set(base, row.sample);
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
    const expectedFrames = row.bases.split(/\s+/).filter(Boolean);
    const matches = appByTruth.get(row.sample) ?? [];
    const problems = [];
    if (matches.length === 0) problems.push("sample missing from the app");
    if (matches.length > 1) problems.push(`split into ${matches.length} app samples`);
    const frames = matches.flatMap((s) => s.frames);
    const seen = frames.map((f) => f.filename);
    const missing = expectedFrames.filter((b) => !seen.includes(b));
    const kept = frames.filter(isKept);
    if (matches.length && missing.length) problems.push(`missing frames: ${missing.join(" ")}`);
    const autoCalls = frames.map((f) => ({ frame: f.filename, call: normalizeCall(f.autoPhases) }));
    const wrongAuto = autoCalls.filter((c) => c.call !== expectedAuto);
    if (wrongAuto.length) {
      problems.push(`auto call ${wrongAuto.map((c) => `${c.frame}=${c.call}`).join(" ")}, expected ${expectedAuto}`);
    }
    const knownGap = expectedAuto !== truth;
    if (strict && knownGap) problems.push(`auto call ${expectedAuto} misses the generated truth ${truth}`);
    const rep = matches.length === 1 ? representativeFrame(matches[0].frames) : null;
    const assigned = rep ? assignedCall(rep.assignment) : null;
    const assignedFrames = frames
      .map((f) => ({ frame: f.filename, call: assignedCall(f.assignment) }))
      .filter((f) => f.call !== null);
    return {
      sample: row.sample,
      sample_id: matches.length === 1 ? matches[0].id : null,
      ok: problems.length === 0,
      problems,
      truth,
      expected_auto: expectedAuto,
      known_gap: knownGap ? row.note || `auto ${expectedAuto} vs truth ${truth}` : null,
      frames: { expected: expectedFrames.length, found: frames.length, kept: kept.length, missing },
      auto_calls: autoCalls,
      representative: rep ? rep.filename : null,
      assigned,
      assigned_matches_truth: assigned === null ? null : assigned === truth,
      assigned_frames: assignedFrames,
    };
  });

  const failed = samples.filter((s) => !s.ok);
  return {
    ok: failed.length === 0 && strays.length === 0,
    gate: strict ? "truth" : "auto",
    summary: {
      samples: samples.length,
      passed: samples.length - failed.length,
      failed: failed.map((s) => s.sample),
      known_gaps: samples.filter((s) => s.known_gap).map((s) => s.sample),
      frames_expected: samples.reduce((n, s) => n + s.frames.expected, 0),
      frames_found: samples.reduce((n, s) => n + s.frames.found, 0),
      assigned: samples.filter((s) => s.assigned !== null).length,
      assigned_off_representative: samples
        .filter((s) => s.assigned === null && s.assigned_frames.length > 0)
        .map((s) => s.sample),
    },
    strays,
    samples,
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
      return { id: s.id, name: s.name, frames };
    }),
  );
}

export async function checkTruth(cfg, flags) {
  const truth = readTruth(cfg, flags.truth);
  if (!truth) throw new UsageError(`no TRUTH.csv at ${flags.truth ?? cfg.demoDir}`, "run `make demo-data` or pass --truth <file>");
  let experimentId = flags.experiment ? Number(flags.experiment) : null;
  if (experimentId === null) {
    const exp = await findDemoExperiment(cfg).catch(() => {
      throw new UsageError(`backend not answering at ${cfg.apiBase}`, "run `bin/verify seed`");
    });
    if (!exp) throw new UsageError("the demo experiment is not ingested", "run `bin/verify seed`");
    if (exp.ingest_status !== "complete") throw new UsageError(`demo experiment ${exp.id} is ${exp.ingest_status}`, "wait for the scan to finish; `bin/verify doctor` shows progress");
    experimentId = exp.id;
  }
  const app = await collectApp(cfg, experimentId);
  return { experiment_id: experimentId, ...compareTruth(truth, app, { strict: Boolean(flags.strict) }) };
}
