import { test } from "node:test";
import assert from "node:assert/strict";
import { readAutosamplerViews, readSample } from "../autosampler_ui.mjs";

const NAME = "wat 01A (B0101)";

function fakePage(stateAt) {
  const clock = { t: 0 };
  const state = () => stateAt(clock.t);
  const element = (text, attributes = {}) => ({
    textContent: async () => text,
    getAttribute: async (name) => attributes[name] ?? null,
  });
  const locator = (selector) => {
    const s = state();
    const found = (n, text = "", attributes = {}) => ({ count: async () => n, first: () => element(text, attributes) });
    if (selector.includes("background-name")) return found(s.background === null ? 0 : 1, s.background ?? "");
    if (selector.includes("background-note")) return found(s.note === null ? 0 : 1, s.note ?? "");
    if (selector.includes("subtraction-toggle")) {
      const visible = s.toggle && !(s.toggle.reserved && selector.includes(":not([aria-hidden])"));
      return found(visible ? 1 : 0, "", { "aria-label": s.toggle?.label ?? "" });
    }
    if (selector.includes("difference-line")) return found(s.difference);
    if (selector.includes("trace-line")) return found(s.trace);
    if (selector.includes("trace-empty")) return found(s.empty === null ? 0 : 1, s.empty ?? "");
    if (selector.includes("workspace-candidates")) return found(s.candidates ?? 0);
    if (selector.includes("peak-key")) return found(s.peaks ?? 0);
    if (selector.includes("custom-index")) return found(s.custom ?? 0);
    throw new Error(`the fake page has no answer for ${selector}`);
  };
  return { locator, waitForTimeout: async (ms) => { clock.t += ms; }, now: () => clock.t, closed: false, close: async function close() { this.closed = true; } };
}

const idle = { background: null, note: null, toggle: null, difference: 0, trace: 0, empty: null };
const options = (page) => ({ now: page.now, timeoutMs: 5000 });

test("a sample whose background trace arrives late is read once Minus is drawn, not before", async () => {
  const page = fakePage((t) => ({
    ...idle,
    background: t < 100 ? "Loading" : `${NAME}run's buffer`,
    note: t < 100 ? null : t < 900 ? "loading" : "",
    toggle: t < 900 ? { reserved: true, label: "Minus " } : { label: `Minus ${NAME}` },
    difference: t < 900 ? 0 : 1,
  }));
  const ui = await readSample(page, { minus: true }, options(page));
  assert.deepEqual([ui.difference_line, ui.minus_pressed, ui.plot_drawn], [1, true, true]);
  assert.equal(ui.background, `${NAME}run's buffer`);
});

test("the reserved toggle shown while the background loads is not a pressed Minus", async () => {
  const page = fakePage(() => ({ ...idle, background: `${NAME}run's buffer`, note: "loading", toggle: { reserved: true, label: "Minus " } }));
  const ui = await readSample(page, { minus: true }, options(page));
  assert.deepEqual([ui.difference_line, ui.minus_pressed, ui.plot_drawn], [0, false, false]);
});

test("a sample that never draws is reported as not drawn after the timeout", async () => {
  const page = fakePage(() => ({ ...idle, background: `${NAME}run's buffer`, note: "no trace" }));
  const ui = await readSample(page, { minus: true }, options(page));
  assert.equal(ui.plot_drawn, false);
  assert.equal(ui.difference_line, 0);
  assert.ok(page.now() >= 5000);
});

test("modules that load with the plot are read after it is drawn", async () => {
  const page = fakePage((t) => ({ ...idle, background: "This is a buffer", trace: t < 1000 ? 0 : 1, candidates: t < 1000 ? 0 : 3 }));
  const ui = await readSample(page, { minus: false }, options(page));
  assert.equal(ui.plot_drawn, true);
  assert.equal(ui.candidates, 3);
});

test("a buffer with no difference line is read once its own trace is drawn", async () => {
  const page = fakePage((t) => ({ ...idle, background: "This is a buffer", trace: t < 300 ? 0 : 1, difference: 0 }));
  const ui = await readSample(page, { minus: false }, options(page));
  assert.deepEqual([ui.plot_drawn, ui.difference_line, ui.minus_pressed], [true, 0, false]);
});

test("a shot with no trace counts as drawn once the plot says so", async () => {
  const page = fakePage((t) => ({ ...idle, background: "No buffer in this run", empty: t < 400 ? "" : "No trace for this image" }));
  const ui = await readSample(page, { minus: false, traceless: true }, options(page));
  assert.equal(ui.plot_drawn, true);
});

function session(page, { gotoFails = false } = {}) {
  const calls = { browserClosed: 0, visited: [] };
  const deps = {
    reachable: async () => ({ ok: true }),
    ensureBrowser: async () => undefined,
    connect: async () => ({ context: { newPage: async () => page }, browser: { close: async () => { calls.browserClosed++; } } }),
    prepareUser: async () => "kept",
    settledGoto: async (_page, _cfg, route) => {
      calls.visited.push(route);
      if (gotoFails) throw new Error("did not settle");
    },
    readSample: async () => ({ background: "This is a buffer", difference_line: 0, minus_pressed: false, plot_drawn: true, candidates: 0, peak_glyphs: 0, custom_index: 0 }),
  };
  return { calls, deps };
}

const cfg = { uiBase: "http://127.0.0.1:1" };

test("each sample is read in a fresh page that is closed afterwards", async () => {
  const page = fakePage(() => idle);
  const { calls, deps } = session(page);
  const out = await readAutosamplerViews(cfg, 1, [{ id: 5, name: "a" }, { id: 6, name: "b" }], deps);
  assert.deepEqual([...out.keys()], [5, 6]);
  assert.deepEqual(calls.visited, ["/experiments/1/workspace/5", "/experiments/1/workspace/6"]);
  assert.equal(page.closed, true);
  assert.equal(calls.browserClosed, 1);
});

test("the fresh page and the connection are closed when a sample fails to load", async () => {
  const page = fakePage(() => idle);
  const { calls, deps } = session(page, { gotoFails: true });
  await assert.rejects(readAutosamplerViews(cfg, 1, [{ id: 5, name: "a" }], deps), /did not settle/);
  assert.equal(page.closed, true);
  assert.equal(calls.browserClosed, 1);
});
