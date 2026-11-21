import { test } from "node:test";
import assert from "node:assert/strict";
import { parseArgs, UsageError } from "../args.mjs";

const spec = { booleans: ["exact", "dry-run"], strings: ["role", "name"] };

test("parses booleans, strings, = form and positionals", () => {
  const r = parseArgs(["GET", "--role", "button", "--name=Save as", "--exact", "/api/x"], spec);
  assert.deepEqual(r.positional, ["GET", "/api/x"]);
  assert.deepEqual(r.flags, { role: "button", name: "Save as", exact: true });
});

test("everything after -- is positional", () => {
  assert.deepEqual(parseArgs(["--", "--role"], spec).positional, ["--role"]);
});

test("unknown option names the options that exist", () => {
  assert.throws(() => parseArgs(["--rol", "x"], spec), (e) => e instanceof UsageError && /--role/.test(e.hint));
});

test("a string option without a value is an error", () => {
  assert.throws(() => parseArgs(["--name"], spec), /needs a value/);
});

test("a boolean option with a value is an error", () => {
  assert.throws(() => parseArgs(["--exact=yes"], spec), /takes no value/);
});

test("parseThrottle reads link speed and latency, and refuses anything else", async () => {
  const { parseThrottle } = await import("../throttle.mjs");
  assert.equal(parseThrottle(undefined), null);
  assert.equal(parseThrottle(""), null);
  assert.deepEqual(parseThrottle("20mbit,30ms"), { mbit: 20, latencyMs: 30 });
  assert.deepEqual(parseThrottle(" 2.5mbit , 0ms "), { mbit: 2.5, latencyMs: 0 });
  assert.throws(() => parseThrottle("fast"), /VERIFY_THROTTLE/);
  assert.throws(() => parseThrottle("0mbit,30ms"), /VERIFY_THROTTLE/);
});
