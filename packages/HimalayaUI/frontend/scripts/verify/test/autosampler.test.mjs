import { test } from "node:test";
import assert from "node:assert/strict";
import { compareAutosampler } from "../autosampler.mjs";

const BUFFER = "lysA 01A (B0101)";
const SAMPLE = "lysA 01B (S0102)";
const ROWS = [
  { sample: BUFFER, stem: "lysA_01A_B0101", kind: "buffer", default_background: "", frames_total: "10", frames_rejected: "0" },
  { sample: SAMPLE, stem: "lysA_01B_S0102", kind: "sample", default_background: BUFFER, frames_total: "10", frames_rejected: "4" },
  { sample: "agbe_S0100", stem: "agbe_S0100", kind: "calibration", default_background: "", frames_total: "", frames_rejected: "" },
];

const CLEAN = { peaks: 0, indices: 0, call_phases: [], members: 0, error: null };

function sample(row, over = {}) {
  const buffer = row.kind === "buffer";
  return {
    id: 1, name: row.sample, mode: "autosampler", data_dir: "/demo/autosampler/data", files: [row.stem],
    frames: row.frames_total === "" ? { total: null, rejected: null } : { total: Number(row.frames_total), rejected: Number(row.frames_rejected) },
    background: { name: row.default_background || null, source: "default", is_buffer: buffer },
    analysis: { ...CLEAN },
    ...over,
  };
}

const good = () => ROWS.map((r) => sample(r));
const options = { dataDir: "/demo/autosampler/data" };
const problemsOf = (r, name) => r.problems.find((p) => p.sample === name)?.problems ?? [];

test("every expected acquisition ingested right, with its background, frames and no analysis, passes", () => {
  const r = compareAutosampler(ROWS, good(), options);
  assert.equal(r.ok, true, JSON.stringify(r.problems));
  assert.deepEqual([r.expected, r.found], [3, 3]);
  assert.deepEqual(r.problems, []);
  assert.deepEqual(r.extras, []);
});

test("a missing, split, multi-image, mislabelled or wrong-mode sample is a problem", () => {
  const missing = compareAutosampler(ROWS, good().slice(1), options);
  assert.equal(missing.ok, false);
  assert.match(problemsOf(missing, BUFFER)[0], /missing/);

  const split = compareAutosampler(ROWS, [...good(), sample(ROWS[0], { id: 9 })], options);
  assert.match(problemsOf(split, BUFFER)[0], /2 samples/);

  const many = good();
  many[1].files = [ROWS[1].stem, ROWS[1].stem + "_0_002"];
  assert.match(problemsOf(compareAutosampler(ROWS, many, options), SAMPLE)[0], /2 images/);

  const other = good();
  other[1].files = ["lysA_01C_S0103"];
  assert.match(problemsOf(compareAutosampler(ROWS, other, options), SAMPLE)[0], /holds lysA_01C_S0103/);

  const capillary = good();
  capillary[0].mode = "capillary";
  assert.match(problemsOf(compareAutosampler(ROWS, capillary, options), BUFFER)[0], /capillary/);

  const elsewhere = good();
  elsewhere[0].data_dir = "/demo/capillary/data";
  assert.match(problemsOf(compareAutosampler(ROWS, elsewhere, options), BUFFER)[0], /\/demo\/capillary\/data/);
});

test("an autosampler sample TRUTH does not list is an extra and fails", () => {
  const r = compareAutosampler(ROWS, [...good(), sample({ ...ROWS[0], sample: "lysA 01Z (S0999)", stem: "lysA_01Z_S0999" }, { id: 5 })], options);
  assert.equal(r.ok, false);
  assert.deepEqual(r.extras, ["lysA 01Z (S0999)"]);
});

test("a sample's default background is its own run's first buffer, and a buffer has none", () => {
  const wrongRun = good();
  wrongRun[1].background.name = "micB 02A (B0108)";
  assert.match(problemsOf(compareAutosampler(ROWS, wrongRun, options), SAMPLE)[0], /background is micB 02A \(B0108\), expected lysA 01A \(B0101\)/);

  const none = good();
  none[1].background.name = null;
  assert.match(problemsOf(compareAutosampler(ROWS, none, options), SAMPLE)[0], /background is none, expected lysA 01A/);

  const chosen = good();
  chosen[1].background.source = "user";
  assert.match(problemsOf(compareAutosampler(ROWS, chosen, options), SAMPLE)[0], /background source is user/);

  const bufferWithOne = good();
  bufferWithOne[0].background.name = BUFFER;
  assert.match(problemsOf(compareAutosampler(ROWS, bufferWithOne, options), BUFFER)[0], /background is lysA 01A \(B0101\), expected none/);

  const flag = good();
  flag[0].background.is_buffer = false;
  assert.match(problemsOf(compareAutosampler(ROWS, flag, options), BUFFER)[0], /not marked as a buffer/);

  const noEntry = good();
  noEntry[1].background = null;
  assert.match(problemsOf(compareAutosampler(ROWS, noEntry, options), SAMPLE)[0], /no entry/);
});

test("frames kept and total match what the generator wrote, and are absent for a calibration shot", () => {
  const rejected = good();
  rejected[1].frames.rejected = 3;
  assert.match(problemsOf(compareAutosampler(ROWS, rejected, options), SAMPLE)[0], /frames 10 total and 3 rejected, expected 10 and 4/);

  const unread = good();
  unread[0].frames = { total: null, rejected: null };
  assert.match(problemsOf(compareAutosampler(ROWS, unread, options), BUFFER)[0], /frames none total and none rejected, expected 10 and 0/);

  const calibration = good();
  calibration[2].frames = { total: 10, rejected: 0 };
  assert.match(problemsOf(compareAutosampler(ROWS, calibration, options), "agbe_S0100")[0], /expected none and none/);
});

test("no automatic analysis ran: no peaks, no indices, no call, no analysis error", () => {
  for (const [field, value, pattern] of [
    ["peaks", 2, /2 peaks/],
    ["indices", 1, /1 indices/],
    ["call_phases", ["Pn3m"], /a call of Pn3m/],
    ["members", 1, /1 ticked indices/],
    ["error", "bad header", /analysis error bad header/],
  ]) {
    const analyzed = good();
    analyzed[1].analysis = { ...CLEAN, [field]: value };
    assert.match(problemsOf(compareAutosampler(ROWS, analyzed, options), SAMPLE)[0], pattern, field);
  }
});

const UI_SAMPLE = { plot_drawn: true, background: `${BUFFER}run's buffer`, difference_line: 1, minus_pressed: true, candidates: 0, peak_glyphs: 0, custom_index: 0 };
const UI_BUFFER = { plot_drawn: true, background: "This is a buffer", difference_line: 0, minus_pressed: false, candidates: 0, peak_glyphs: 0, custom_index: 0 };
const UI_CALIBRATION = { plot_drawn: true, background: "No buffer in this run", difference_line: 0, minus_pressed: false, candidates: 0, peak_glyphs: 0, custom_index: 0 };

function withUi() {
  const list = good();
  list[0].ui = { ...UI_BUFFER };
  list[1].ui = { ...UI_SAMPLE };
  list[2].ui = { ...UI_CALIBRATION };
  return list;
}

test("with the workspace read, a sample opens in Minus on its run's buffer and a buffer says so", () => {
  assert.equal(compareAutosampler(ROWS, withUi(), options).ok, true);

  const raw = withUi();
  raw[1].ui = { ...UI_SAMPLE, difference_line: 0, minus_pressed: false };
  const problems = problemsOf(compareAutosampler(ROWS, raw, options), SAMPLE);
  assert.ok(problems.some((p) => /not drawn in Minus/.test(p)), problems.join("; "));

  const chosen = withUi();
  chosen[1].ui = { ...UI_SAMPLE, background: `${BUFFER}chosen` };
  assert.match(problemsOf(compareAutosampler(ROWS, chosen, options), SAMPLE)[0], /Background module reads "lysA 01A \(B0101\)chosen"/);

  const wrongBuffer = withUi();
  wrongBuffer[1].ui = { ...UI_SAMPLE, background: "micB 02A (B0108)run's buffer" };
  assert.match(problemsOf(compareAutosampler(ROWS, wrongBuffer, options), SAMPLE)[0], /Background module reads/);

  const bufferSubtracted = withUi();
  bufferSubtracted[0].ui = { ...UI_BUFFER, background: BUFFER + "run's buffer" };
  assert.match(problemsOf(compareAutosampler(ROWS, bufferSubtracted, options), BUFFER)[0], /Background module reads/);
});

test("with the workspace read, the peak, call and indexing modules are absent", () => {
  for (const [field, pattern] of [["candidates", /call candidates/], ["peak_glyphs", /peak markers/], ["custom_index", /custom index panel/]]) {
    const list = withUi();
    list[1].ui = { ...UI_SAMPLE, [field]: 1 };
    assert.match(problemsOf(compareAutosampler(ROWS, list, options), SAMPLE)[0], pattern, field);
  }
});

test("a plot that never drew is a problem, whatever the absence checks read", () => {
  const list = withUi();
  list[1].ui = { ...UI_SAMPLE, plot_drawn: false, difference_line: 0, minus_pressed: false };
  const problems = problemsOf(compareAutosampler(ROWS, list, options), SAMPLE);
  assert.match(problems[0], /plot was not drawn/);
  assert.equal(problems.some((p) => /Minus/.test(p)), false);
});

test("a truth file without the frame count or background columns fails instead of matching a null", () => {
  const stripped = ROWS.map(({ frames_total: _t, frames_rejected: _r, ...row }) => row);
  const r = compareAutosampler(stripped, good(), options);
  assert.equal(r.ok, false);
  const columns = r.problems.find((p) => p.sample === "TRUTH_autosampler.csv");
  assert.match(columns.problems[0], /missing the columns frames_total, frames_rejected/);

  const noBackground = ROWS.map(({ default_background: _b, ...row }) => row);
  assert.match(compareAutosampler(noBackground, good(), options).problems.find((p) => p.sample === "TRUTH_autosampler.csv").problems[0], /default_background/);
});

test("an empty frame count is legal only for the calibration shot", () => {
  const rows = ROWS.map((r) => (r.sample === SAMPLE ? { ...r, frames_total: "", frames_rejected: "" } : r));
  const list = rows.map((r) => sample(r));
  assert.match(problemsOf(compareAutosampler(rows, list, options), SAMPLE)[0], /no frame counts for this sample/);
  assert.equal(compareAutosampler(ROWS, good(), options).ok, true);
});
