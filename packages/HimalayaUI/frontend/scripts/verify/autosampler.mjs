import { getJson } from "./http.mjs";

export const AUTOSAMPLER_MODE = "autosampler";

export function setupMode(setups, sample) {
  return setups.find((s) => s.id === sample.setup_id)?.mode ?? null;
}

async function analysisOf(api, sample, exposure) {
  const [peaks, indices, call] = await Promise.all([
    exposure ? getJson(api, `/api/exposures/${exposure.id}/peaks`) : [],
    exposure ? getJson(api, `/api/exposures/${exposure.id}/indices`) : [],
    getJson(api, `/api/samples/${sample.id}/call`),
  ]);
  return {
    peaks: peaks.length,
    indices: indices.length,
    call_phases: [...(call.call_phases ?? [])],
    members: (call.members ?? []).length,
    error: exposure?.analysis_error ?? null,
  };
}

export async function collectAutosampler(api, experimentId, samples, setups, rows) {
  const expected = new Set(rows.map((r) => r.sample));
  const wanted = samples.filter((s) => setupMode(setups, s) === AUTOSAMPLER_MODE || expected.has(s.name));
  const compares = await getJson(api, `/api/experiments/${experimentId}/compares`);
  const backgrounds = new Map((compares.backgrounds ?? []).map((b) => [b.sample_id, b]));
  const names = new Map(samples.map((s) => [s.id, s.name]));
  return Promise.all(wanted.map(async (s) => {
    const setup = setups.find((x) => x.id === s.setup_id);
    const exposures = (await getJson(api, `/api/samples/${s.id}/exposures`)).filter((e) => e.kind !== "combined");
    const only = exposures.length === 1 ? exposures[0] : null;
    const background = backgrounds.get(s.id);
    return {
      id: s.id,
      name: s.name,
      mode: setup?.mode ?? null,
      data_dir: setup?.data_dir ?? null,
      files: exposures.map((e) => e.filename),
      exposure_id: only?.id ?? null,
      frames: only ? { total: only.frames_total ?? null, rejected: only.frames_rejected ?? null } : null,
      background: background
        ? { name: background.ref_sample_id == null ? null : names.get(background.ref_sample_id) ?? `sample ${background.ref_sample_id}`, source: background.source, is_buffer: background.is_buffer }
        : null,
      analysis: await analysisOf(api, s, only),
    };
  }));
}

const orNone = (v) => (v === null || v === undefined || v === "" ? "none" : v);
const count = (v) => (v === "" || v === undefined ? null : Number(v));

function backgroundProblems(row, s) {
  if (!s.background) return ["no entry in the compares' backgrounds"];
  const expected = row.default_background || null;
  const out = [];
  if (s.background.name !== expected) out.push(`background is ${orNone(s.background.name)}, expected ${orNone(expected)}`);
  if (s.background.source !== "default") out.push(`background source is ${s.background.source}, expected default`);
  const buffer = row.kind === "buffer";
  if (s.background.is_buffer !== buffer) out.push(buffer ? "not marked as a buffer" : "marked as a buffer");
  return out;
}

function frameProblems(row, s) {
  if (!s.frames) return [];
  const total = count(row.frames_total);
  const rejected = count(row.frames_rejected);
  if (s.frames.total === total && s.frames.rejected === rejected) return [];
  return [`frames ${orNone(s.frames.total)} total and ${orNone(s.frames.rejected)} rejected, expected ${orNone(total)} and ${orNone(rejected)}`];
}

function analysisProblems(s) {
  const a = s.analysis;
  const out = [];
  if (a.peaks) out.push(`has ${a.peaks} peaks, but no analysis runs here`);
  if (a.indices) out.push(`has ${a.indices} indices, but no analysis runs here`);
  if (a.call_phases.length) out.push(`has a call of ${a.call_phases.join("+")}`);
  if (a.members) out.push(`has ${a.members} ticked indices`);
  if (a.error) out.push(`has an analysis error ${a.error}`);
  return out;
}

function expectedBackgroundText(row) {
  if (row.kind === "buffer") return "This is a buffer";
  return row.default_background ? `${row.default_background}run's buffer` : "No buffer in this run";
}

function uiProblems(row, ui) {
  const out = [];
  const expected = expectedBackgroundText(row);
  if (ui.background !== expected) out.push(`Background module reads "${ui.background}", expected "${expected}"`);
  const minus = Boolean(row.default_background) && row.kind !== "buffer";
  if (minus && !(ui.difference_line === 1 && ui.minus_pressed)) out.push(`not drawn in Minus: ${ui.difference_line} difference lines, Minus ${ui.minus_pressed ? "" : "not "}pressed`);
  if (!minus && ui.difference_line) out.push("draws a difference line without a background");
  if (ui.candidates) out.push("shows the call candidates");
  if (ui.peak_glyphs) out.push("shows peak markers");
  if (ui.custom_index) out.push("shows the custom index panel");
  return out;
}

function rowProblems(row, matches, dataDir) {
  if (matches.length === 0) return ["sample missing from the app"];
  if (matches.length > 1) return [`${matches.length} samples carry this name`];
  const [s] = matches;
  const problems = [];
  if (s.mode !== AUTOSAMPLER_MODE) problems.push(`sits in a ${s.mode ?? "mode-less"} setup, not an autosampler one`);
  else if (dataDir && s.data_dir !== dataDir) problems.push(`sits in the setup ${s.data_dir}, expected ${dataDir}`);
  if (s.files.length !== 1) problems.push(`has ${s.files.length} images, expected one`);
  else if (s.files[0] !== row.stem) problems.push(`holds ${s.files[0]}, expected ${row.stem}`);
  if (s.mode === AUTOSAMPLER_MODE) {
    if (s.ui) problems.push(...uiProblems(row, s.ui));
    problems.push(...backgroundProblems(row, s), ...frameProblems(row, s), ...analysisProblems(s));
  }
  return problems;
}

export function compareAutosampler(rows, samples, { dataDir = null } = {}) {
  const byName = new Map();
  for (const s of samples) byName.set(s.name, [...(byName.get(s.name) ?? []), s]);
  const problems = rows
    .map((row) => ({ sample: row.sample, problems: rowProblems(row, byName.get(row.sample) ?? [], dataDir) }))
    .filter((p) => p.problems.length);
  const known = new Set(rows.map((r) => r.sample));
  const extras = samples.filter((s) => s.mode === AUTOSAMPLER_MODE && !known.has(s.name)).map((s) => s.name);
  return {
    ok: problems.length === 0 && extras.length === 0,
    gate: "ingestion, default background, frame counts, no analysis, and with --ui the workspace's Background module and Minus view",
    expected: rows.length,
    found: rows.filter((r) => (byName.get(r.sample) ?? []).length === 1).length,
    problems,
    extras,
  };
}
