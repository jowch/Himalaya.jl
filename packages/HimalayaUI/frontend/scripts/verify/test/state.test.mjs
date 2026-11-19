import { test } from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import { FRONTEND } from "../config.mjs";
import { STATE_KEY, STATE_VERSION } from "../page.mjs";

test("the seeded localStorage key and version match the app's persisted store", () => {
  const src = fs.readFileSync(path.join(FRONTEND, "src/state.ts"), "utf8");
  assert.match(src, new RegExp(`LS_KEY = "${STATE_KEY}"`));
  const persist = src.slice(src.indexOf("name: LS_KEY"));
  assert.match(persist, new RegExp(`version: ${STATE_VERSION},`));
  assert.match(src, /tutorialSeen: boolean/);
  assert.match(src, /username: string \| undefined/);
});
