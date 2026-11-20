import fs from "node:fs";
import { loadConfig } from "../frontend/scripts/verify/config.mjs";
import { getJson, request } from "../frontend/scripts/verify/http.mjs";

const read = (file) => JSON.parse(fs.readFileSync(file, "utf8"));
const report = (out) => {
  console.log(JSON.stringify(out, null, 2));
  process.exit(out.ok ? 0 : 1);
};
const sampleOf = (result, label) => {
  const s = (result.samples ?? []).find((x) => x.sample === label);
  if (!s) report({ ok: false, error: `no sample ${label} in the check-truth output` });
  return s;
};

async function shownExposure(api, s, label) {
  const exposures = await getJson(api, `/api/samples/${s.sample_id}/exposures`);
  const exposure = exposures.find((e) => e.filename === s.shown.frame);
  if (!exposure) report({ ok: false, error: `${label}'s shown frame ${s.shown.frame} is not among its exposures` });
  return exposure;
}

async function confirm(file, label) {
  const s = sampleOf(read(file), label);
  const call = s.expected_auto;
  if (call === "none" || call.includes("+")) report({ ok: false, error: `${label}'s auto call is ${call}; pick a single-phase sample` });
  const api = loadConfig().apiBase;
  const exposure = await shownExposure(api, s, label);
  const indices = await getJson(api, `/api/exposures/${exposure.id}/indices`);
  const index = indices.filter((i) => i.phase === call).sort((a, b) => b.score - a.score)[0];
  if (!index) report({ ok: false, error: `${s.shown.frame} has no ${call} index to confirm` });
  const res = await request(api, "POST", `/api/exposures/${exposure.id}/assignment/members`, { data: { index_id: index.id } });
  report({ ok: res.ok, rule: "a person confirms the auto call on one sample", sample: label, frame: s.shown.frame, exposure_id: exposure.id, index_id: index.id, call, status: res.status });
}

function confirmed(file, label) {
  const s = sampleOf(read(file), label);
  const want = s.expected_auto;
  const assigned = s.assigned_frames.map((f) => f.call);
  const ok = s.shown?.status === "phase" && s.shown.call === want && assigned.includes(want);
  report({ ok, rule: `${label} shows the confirmed ${want} call, in the API and in the corpus`, shown: s.shown, assigned_frames: s.assigned_frames });
}

async function target(file, label) {
  const s = sampleOf(read(file), label);
  const exposure = await shownExposure(loadConfig().apiBase, s, label);
  console.log(`${s.sample_id} ${exposure.id}`);
}

async function custom(file, label, phase, want) {
  const s = sampleOf(read(file), label);
  const api = loadConfig().apiBase;
  const exposure = await shownExposure(api, s, label);
  const indices = await getJson(api, `/api/exposures/${exposure.id}/indices`);
  const customs = indices.filter((i) => i.kind === "speculative" && i.phase === phase);
  const assignment = await getJson(api, `/api/exposures/${exposure.id}/assignment`);
  const members = assignment?.members ?? [];
  const assigned = customs.filter((i) => members.includes(i.id)).map((i) => i.id);
  const ok = customs.length === 1 && (want === "on" ? assigned.length === 1 : assigned.length === 0);
  report({ ok, rule: `a custom ${phase} index added in the workspace is stored once and is ${want === "on" ? "in" : "out of"} the frame's call`, sample: label, exposure_id: exposure.id, custom_ids: customs.map((i) => i.id), members });
}

function unchanged(beforeFile, afterFile) {
  const pick = (r) => Object.fromEntries((r.samples ?? []).map((s) => [s.sample, { shown: s.shown, views: s.views, auto_calls: s.auto_calls, assigned_frames: s.assigned_frames }]));
  const before = pick(read(beforeFile));
  const after = pick(read(afterFile));
  const changed = Object.keys({ ...before, ...after })
    .filter((k) => JSON.stringify(before[k]) !== JSON.stringify(after[k]))
    .map((k) => ({ sample: k, before: before[k] ?? null, after: after[k] ?? null }));
  report({ ok: Object.keys(before).length > 0 && changed.length === 0, rule: "re-indexing from scratch changes nothing a person sees, and keeps their confirmed call", samples: Object.keys(before).length, changed });
}

const peakName = (p) => `${p.source === "manual" ? "Manual" : "Auto"} peak at q = ${p.q.toFixed(4)}${p.excluded ? " (excluded)" : ""}`;

function gapQ(peaks) {
  const qs = peaks.map((p) => p.q).sort((a, b) => a - b);
  let best = null;
  for (let i = 1; i < qs.length; i++) {
    if (!best || qs[i] - qs[i - 1] > best[1] - best[0]) best = [qs[i - 1], qs[i]];
  }
  return best ? (best[0] + best[1]) / 2 : null;
}

async function addManual(api, exposureId, peaks) {
  const q = gapQ(peaks);
  if (q === null) report({ ok: false, error: `exposure ${exposureId} has fewer than two peaks to place a manual peak between` });
  const res = await request(api, "POST", `/api/exposures/${exposureId}/peaks`, { data: { q } });
  if (!res.ok) report({ ok: false, error: `adding a manual peak to exposure ${exposureId} returned ${res.status}`, body: res.body });
  return res.body;
}

async function twin(file) {
  const api = loadConfig().apiBase;
  const frames = [];
  for (const s of read(file).samples ?? []) {
    if (s.assigned_frames?.length) continue;
    const exposures = await getJson(api, `/api/samples/${s.sample_id}/exposures`);
    for (const [i, exposure] of exposures.entries()) {
      if (exposure.status === "rejected") continue;
      const peaks = await getJson(api, `/api/exposures/${exposure.id}/peaks`);
      if (peaks.length >= 2 && peaks.some((p) => p.source === "auto")) frames.push({ s, exposure, frame: `Frame ${i + 1}`, peaks });
    }
  }
  if (!frames.length) report({ ok: false, error: "no kept frame of an unassigned sample has two or more peaks, one of them auto" });
  const lowest = (f) => Math.min(...f.peaks.filter((p) => p.source === "auto").map((p) => p.id));
  frames.sort((a, b) => lowest(a) - lowest(b));
  const autoWith = (id) => frames.find((f) => f.peaks.some((p) => p.source === "auto" && p.id === id));
  let host = frames[0];
  let manual = await addManual(api, host.exposure.id, host.peaks);
  if (autoWith(manual.id) !== host) {
    const del = await request(api, "DELETE", `/api/peaks/manual/${manual.id}`);
    if (!del.ok) report({ ok: false, error: `removing the probe manual peak ${manual.id} returned ${del.status}` });
    host = autoWith(manual.id + 1);
    if (!host) report({ ok: false, error: `no kept frame has an auto peak with id ${manual.id + 1}, the next manual peak id, so the check could not make the ids collide and would pass vacuously` });
    manual = await addManual(api, host.exposure.id, host.peaks);
  }
  const auto = host.peaks.find((p) => p.source === "auto" && p.id === manual.id);
  if (!auto) report({ ok: false, error: `manual peak ${manual.id} on exposure ${host.exposure.id} shares no id with an auto peak there`, manual });
  console.log(JSON.stringify({ sample_id: host.s.sample_id, exposure_id: host.exposure.id, frame: host.frame, id: manual.id, manual: peakName(manual), auto: peakName(auto), excluded: peakName({ ...auto, excluded: true }) }));
}

async function twinState(exposureId, id, want) {
  const peaks = await getJson(loadConfig().apiBase, `/api/exposures/${exposureId}/peaks`);
  const same = peaks.filter((p) => p.id === Number(id));
  const manual = same.filter((p) => p.source === "manual");
  const auto = same.filter((p) => p.source === "auto");
  const expected = { removed: false, excluded: true, restored: false }[want];
  const ok = manual.length === 0 && auto.length === 1 && auto[0].excluded === expected;
  report({ ok, rule: `clicking a peak acts on that peak alone, when an auto and a manual peak share id ${id} (${want})`, exposure_id: Number(exposureId), with_id: same });
}

const [cmd, ...args] = process.argv.slice(2);
const commands = { confirm, confirmed, unchanged, target, custom, twin, "twin-state": twinState };
if (!commands[cmd]) report({ ok: false, error: `usage: verify_truth.mjs ${Object.keys(commands).join("|")} ...` });
await commands[cmd](...args);
