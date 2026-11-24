import { getJson } from "./http.mjs";

export const AUTOSAMPLER_MODE = "autosampler";

export function setupMode(setups, sample) {
  return setups.find((s) => s.id === sample.setup_id)?.mode ?? null;
}

export async function collectAutosampler(api, samples, setups, rows) {
  const expected = new Set(rows.map((r) => r.sample));
  const wanted = samples.filter((s) => setupMode(setups, s) === AUTOSAMPLER_MODE || expected.has(s.name));
  return Promise.all(wanted.map(async (s) => {
    const setup = setups.find((x) => x.id === s.setup_id);
    const exposures = await getJson(api, `/api/samples/${s.id}/exposures`);
    return {
      id: s.id,
      name: s.name,
      mode: setup?.mode ?? null,
      data_dir: setup?.data_dir ?? null,
      files: exposures.filter((e) => e.kind !== "combined").map((e) => e.filename),
    };
  }));
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
    gate: "ingestion only: one sample and one image per acquisition, in the autosampler setup",
    expected: rows.length,
    found: rows.filter((r) => (byName.get(r.sample) ?? []).length === 1).length,
    problems,
    extras,
  };
}
