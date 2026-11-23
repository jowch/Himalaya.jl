#!/usr/bin/env node
import fs from "node:fs";

const [mode, port, file] = process.argv.slice(2);
if (!["review", "check"].includes(mode) || !port || !file) {
  console.error("usage: upgrade_history.mjs review|check <api port> <calls.json>");
  process.exit(2);
}
const base = `http://127.0.0.1:${port}`;

async function api(method, path, body) {
  const res = await fetch(base + path, {
    method,
    headers: { "Content-Type": "application/json", "X-Username": "upgrade-drill" },
    body: body === undefined ? undefined : JSON.stringify(body),
  });
  const text = await res.text();
  if (!res.ok) throw new Error(`${method} ${path} → ${res.status}: ${text.slice(0, 300)}`);
  return text ? JSON.parse(text) : null;
}

async function demoSamples() {
  const experiments = await api("GET", "/api/experiments");
  const demo = experiments.find((e) => e.name === "Demo (synthetic)") ?? experiments[0];
  if (!demo) throw new Error("no experiment to review");
  const samples = await api("GET", `/api/samples?experiment_id=${demo.id}`);
  return samples.sort((a, b) => a.id - b.id);
}

const imagesOf = async (sid) =>
  (await api("GET", `/api/samples/${sid}/exposures`))
    .filter((e) => (e.kind ?? "file") === "file" && e.status !== "rejected")
    .sort((a, b) => a.id - b.id);

async function topCandidate(eid) {
  const indices = await api("GET", `/api/exposures/${eid}/indices`);
  return [...indices].sort((a, b) => (b.score ?? 0) - (a.score ?? 0))[0];
}

async function callOn(eid) {
  const ix = await topCandidate(eid);
  if (!ix) return null;
  await api("POST", `/api/exposures/${eid}/assignment/members`, { index_id: ix.id });
  return ix.phase;
}

async function addPeak(eid) {
  const peaks = await api("GET", `/api/exposures/${eid}/peaks`);
  const qs = peaks.map((p) => p.q).sort((a, b) => a - b);
  if (qs.length < 2) return false;
  const q = qs[qs.length - 1] * 1.07;
  await api("POST", `/api/exposures/${eid}/peaks`, { q });
  return true;
}

const SCENARIOS = [
  ["called on the oldest image, then a peak added on the newest", 2, async (imgs) => (await callOn(imgs[0].id)) && addPeak(imgs.at(-1).id)],
  ["called on the newest image, which is then dropped", 2, async (imgs) => {
    if (!(await callOn(imgs.at(-1).id))) return false;
    await api("PATCH", `/api/exposures/${imgs.at(-1).id}/status`, { status: "rejected" });
    return true;
  }],
  ["called on two images", 2, async (imgs) => Boolean((await callOn(imgs[0].id)) && (await callOn(imgs[1].id)))],
  ["marked no phase", 1, async (imgs, sid) => {
    await api("POST", `/api/samples/${sid}/call/state`, { state: "form_factor" });
    return true;
  }],
];

async function calls() {
  const out = {};
  for (const s of await demoSamples()) {
    const c = await api("GET", `/api/samples/${s.id}/call`);
    out[s.name] = { state: c.state, phases: [...(c.call_phases ?? [])].sort().join("+") || "none" };
  }
  return out;
}

async function callThenSelectOther(samples) {
  for (const s of samples) {
    const imgs = await imagesOf(s.id);
    const phases = new Map();
    for (const e of imgs) phases.set(e.id, new Set((await api("GET", `/api/exposures/${e.id}/indices`)).map((ix) => ix.phase)));
    for (const a of imgs) {
      const ix = await topCandidate(a.id);
      if (!ix) continue;
      const b = imgs.find((e) => e.id !== a.id && !phases.get(e.id).has(ix.phase));
      if (!b) continue;
      await api("POST", `/api/exposures/${a.id}/assignment/members`, { index_id: ix.id });
      await api("PUT", `/api/samples/${s.id}/selection`, { exposure_ids: [b.id] });
      return { sample: s.name, what: `called ${ix.phase} on one image, then selected only an image without it` };
    }
  }
  throw new Error("no sample has an image whose top candidate's phase another kept image lacks, so the drill can't reach a call the selected curve doesn't offer");
}

if (mode === "review") {
  const samples = await demoSamples();
  const first = await callThenSelectOther(samples);
  const done = [first];
  let next = 0;
  for (const s of samples) {
    if (s.name === first.sample) continue;
    if (next >= SCENARIOS.length) break;
    const imgs = await imagesOf(s.id);
    const [what, least, run] = SCENARIOS[next];
    if (imgs.length < least) continue;
    if (await run(imgs, s.id)) {
      done.push({ sample: s.name, what });
      next++;
    }
  }
  if (next < SCENARIOS.length) throw new Error(`only ${next} of ${SCENARIOS.length} review scenarios found a sample with enough kept images and a candidate`);
  const before = await calls();
  fs.writeFileSync(file, JSON.stringify({ reviewed: done, calls: before }, null, 2));
  console.log(JSON.stringify({ reviewed: done }, null, 2));
} else {
  const { calls: before } = JSON.parse(fs.readFileSync(file, "utf8"));
  const after = await calls();
  const changed = Object.keys(before)
    .filter((name) => JSON.stringify(before[name]) !== JSON.stringify(after[name]))
    .map((name) => ({ sample: name, before: before[name], after: after[name] ?? null }));
  if (changed.length) {
    console.log(JSON.stringify({ ok: false, reason: "the upgrade changed these samples' calls", changed }, null, 2));
    process.exit(1);
  }
  console.log(JSON.stringify({ ok: true, samples: Object.keys(before).length }));
}
