import { test } from "node:test";
import assert from "node:assert/strict";
import { requireSampleRows, waitForSampleRows } from "../views.mjs";

function pageShowing(counts) {
  let reads = 0;
  return {
    reads: () => reads,
    locator: () => ({ count: async () => counts[Math.min(reads++, counts.length - 1)] }),
    waitForTimeout: async () => {},
  };
}

test("the workspace is read only once its sample list has every row", async () => {
  const page = pageShowing([0, 0, 4, 13, 13]);
  assert.equal(await waitForSampleRows(page, 13), 13);
  assert.equal(page.reads(), 4);
});

test("a list that never fills is reported as it is, after the timeout", async () => {
  const page = pageShowing([0, 2]);
  assert.equal(await waitForSampleRows(page, 13, { timeoutMs: -1 }), 0);
  assert.equal(await waitForSampleRows(pageShowing([2]), 13, { timeoutMs: -1 }), 2);
});

test("a list still short after the timeout is an error naming how many rows were expected and seen", async () => {
  await assert.rejects(requireSampleRows(pageShowing([2]), 13, { timeoutMs: -1 }), /listed 2 of 13 samples/);
  assert.equal(await requireSampleRows(pageShowing([13]), 13, { timeoutMs: -1 }), 13);
});
