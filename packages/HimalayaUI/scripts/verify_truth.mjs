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

async function confirm(file, label) {
  const s = sampleOf(read(file), label);
  const call = s.expected_auto;
  if (call === "none" || call.includes("+")) report({ ok: false, error: `${label}'s auto call is ${call}; pick a single-phase sample` });
  const api = loadConfig().apiBase;
  const exposures = await getJson(api, `/api/samples/${s.sample_id}/exposures`);
  const exposure = exposures.find((e) => e.filename === s.shown.frame);
  if (!exposure) report({ ok: false, error: `${label}'s shown frame ${s.shown.frame} is not among its exposures` });
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

function unchanged(beforeFile, afterFile) {
  const pick = (r) => Object.fromEntries((r.samples ?? []).map((s) => [s.sample, { shown: s.shown, views: s.views, auto_calls: s.auto_calls, assigned_frames: s.assigned_frames }]));
  const before = pick(read(beforeFile));
  const after = pick(read(afterFile));
  const changed = Object.keys({ ...before, ...after })
    .filter((k) => JSON.stringify(before[k]) !== JSON.stringify(after[k]))
    .map((k) => ({ sample: k, before: before[k] ?? null, after: after[k] ?? null }));
  report({ ok: Object.keys(before).length > 0 && changed.length === 0, rule: "re-indexing from scratch changes nothing a person sees, and keeps their confirmed call", samples: Object.keys(before).length, changed });
}

const [cmd, ...args] = process.argv.slice(2);
const commands = { confirm, confirmed, unchanged };
if (!commands[cmd]) report({ ok: false, error: `usage: verify_truth.mjs ${Object.keys(commands).join("|")} ...` });
await commands[cmd](...args);
