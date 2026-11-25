import { UsageError } from "./args.mjs";
import { ensureBrowser, connect } from "./browser.mjs";
import { reachable } from "./http.mjs";
import { prepareUser } from "./page.mjs";
import { WORKSPACE_ROUTE, settledGoto } from "./views.mjs";

const BACKGROUND = "[data-testid='background-name']";
const NOTE = "[data-testid='background-note']";
const TOGGLE = "[data-testid='workspace-subtraction-toggle']:not([aria-hidden]) button[aria-pressed='true']";
const DIFFERENCE = "[data-role='difference-line']";
const TRACE = "[data-role='trace-line']";
const NO_TRACE = "[data-testid='workspace-trace-empty']";
const SETTLE_MS = 10000;
const POLL_MS = 100;

async function textOf(page, selector) {
  const el = page.locator(selector);
  return (await el.count()) ? ((await el.first().textContent()) ?? "").trim() : null;
}

async function minusPressed(page) {
  const pressed = page.locator(TOGGLE);
  return (await pressed.count()) > 0 && /^Minus/.test((await pressed.first().getAttribute("aria-label")) ?? "");
}

async function plotDrawn(page, expect) {
  if ((await page.locator(DIFFERENCE).count()) || (await page.locator(TRACE).count())) return true;
  return expect.traceless === true && Boolean(await textOf(page, NO_TRACE));
}

async function waitFor(page, ready, { now, timeoutMs }) {
  const until = now() + timeoutMs;
  while (!(await ready()) && now() < until) await page.waitForTimeout(POLL_MS);
  return ready();
}

export async function readSample(page, expect = {}, { now = Date.now, timeoutMs = SETTLE_MS } = {}) {
  const wait = { now, timeoutMs };
  await waitFor(page, async () => {
    const text = await textOf(page, BACKGROUND);
    return text !== null && text !== "Loading";
  }, wait);
  if (expect.minus) {
    await waitFor(page, async () => {
      if ((await page.locator(DIFFERENCE).count()) > 0) return true;
      const note = await textOf(page, NOTE);
      return note !== null && note !== "loading";
    }, wait);
  }
  const drawn = await waitFor(page, () => plotDrawn(page, expect), wait);
  return {
    plot_drawn: drawn,
    background: await textOf(page, BACKGROUND),
    difference_line: await page.locator(DIFFERENCE).count(),
    minus_pressed: await minusPressed(page),
    candidates: await page.locator("[data-testid='workspace-candidates']").count(),
    peak_glyphs: await page.locator("[data-role='peak-glyph'] [data-peak-key]").count(),
    custom_index: await page.locator("[data-testid='workspace-custom-index']").count(),
  };
}

const DEFAULT_DEPS = { reachable, ensureBrowser, connect, prepareUser, settledGoto, readSample };

export async function readAutosamplerViews(cfg, experimentId, samples, deps = DEFAULT_DEPS) {
  const ui = await deps.reachable(cfg.uiBase, "/");
  if (!ui?.ok) {
    throw new UsageError(`nothing is serving ${cfg.uiBase}`, "--ui reads the workspace in the browser; run `bin/verify seed`");
  }
  await deps.ensureBrowser(cfg);
  const session = await deps.connect(cfg);
  const page = await session.context.newPage();
  const out = new Map();
  try {
    await deps.prepareUser(page, cfg, { username: "verify", onboarding: false });
    const held = new Set(samples.flatMap((s) => s.traceless ?? []));
    for (const s of samples) {
      await deps.settledGoto(page, cfg, WORKSPACE_ROUTE(experimentId, s.id), held);
      out.set(s.id, await deps.readSample(page, { minus: s.minus === true, traceless: (s.traceless ?? []).length > 0 }));
    }
  } finally {
    await page.close().catch(() => undefined);
    await session.browser.close().catch(() => undefined);
  }
  return out;
}
