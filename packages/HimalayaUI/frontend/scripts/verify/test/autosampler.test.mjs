import { test } from "node:test";
import assert from "node:assert/strict";
import { compareAutosampler } from "../autosampler.mjs";

const ROWS = [
  { sample: "lysA 01A (B0101)", stem: "lysA_01A_B0101", kind: "buffer" },
  { sample: "lysA 01B (S0102)", stem: "lysA_01B_S0102", kind: "sample" },
];

const sample = (name, stem, over = {}) => ({ id: 1, name, mode: "autosampler", data_dir: "/demo/autosampler/data", files: [stem], ...over });
const good = () => [sample(ROWS[0].sample, ROWS[0].stem), sample(ROWS[1].sample, ROWS[1].stem)];
const options = { dataDir: "/demo/autosampler/data" };

test("every expected acquisition ingested as one single-image sample passes", () => {
  const r = compareAutosampler(ROWS, good(), options);
  assert.equal(r.ok, true);
  assert.deepEqual([r.expected, r.found], [2, 2]);
  assert.deepEqual(r.problems, []);
  assert.deepEqual(r.extras, []);
});

test("a missing, split, multi-image, mislabelled or wrong-mode sample is a problem", () => {
  const missing = compareAutosampler(ROWS, good().slice(1), options);
  assert.equal(missing.ok, false);
  assert.match(missing.problems[0].problems[0], /missing/);

  const split = compareAutosampler(ROWS, [...good(), sample(ROWS[0].sample, ROWS[0].stem, { id: 9 })], options);
  assert.match(split.problems[0].problems[0], /2 samples/);

  const many = good();
  many[1].files = [ROWS[1].stem, ROWS[1].stem + "_0_002"];
  assert.match(compareAutosampler(ROWS, many, options).problems[0].problems[0], /2 images/);

  const other = good();
  other[1].files = ["lysA_01C_S0103"];
  assert.match(compareAutosampler(ROWS, other, options).problems[0].problems[0], /holds lysA_01C_S0103/);

  const capillary = good();
  capillary[0].mode = "capillary";
  assert.match(compareAutosampler(ROWS, capillary, options).problems[0].problems[0], /capillary/);

  const elsewhere = good();
  elsewhere[0].data_dir = "/demo/capillary/data";
  assert.match(compareAutosampler(ROWS, elsewhere, options).problems[0].problems[0], /\/demo\/capillary\/data/);
});

test("an autosampler sample TRUTH does not list is an extra and fails", () => {
  const r = compareAutosampler(ROWS, [...good(), sample("lysA 01Z (S0999)", "lysA_01Z_S0999", { id: 5 })], options);
  assert.equal(r.ok, false);
  assert.deepEqual(r.extras, ["lysA 01Z (S0999)"]);
});
