import { test } from "node:test";
import assert from "node:assert/strict";
import { parseCsv } from "../csv.mjs";

test("parses quoted fields with commas, quotes and CRLF", () => {
  const rows = parseCsv('sample,note\r\nA,"x, ""y"""\r\nB,plain\n');
  assert.deepEqual(rows, [
    { sample: "A", note: 'x, "y"' },
    { sample: "B", note: "plain" },
  ]);
});

test("missing trailing fields read as empty strings", () => {
  assert.deepEqual(parseCsv("a,b,c\n1,2\n"), [{ a: "1", b: "2", c: "" }]);
});

test("empty input gives no rows", () => {
  assert.deepEqual(parseCsv(""), []);
});
