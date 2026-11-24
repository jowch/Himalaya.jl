import { UsageError } from "./args.mjs";
import { ensureBrowser, connect } from "./browser.mjs";
import { reachable } from "./http.mjs";
import { prepareUser } from "./page.mjs";
import { WORKSPACE_ROUTE, settledGoto } from "./views.mjs";

const BACKGROUND = "[data-testid='background-name']";
const SETTLE_MS = 10000;

async function textOf(page, selector) {
  const el = page.locator(selector);
  return (await el.count()) ? ((await el.first().textContent()) ?? "").trim() : null;
}

async function readSample(page) {
  const until = Date.now() + SETTLE_MS;
  let background = await textOf(page, BACKGROUND);
  while ((background === null || background === "Loading") && Date.now() < until) {
    await page.waitForTimeout(100);
    background = await textOf(page, BACKGROUND);
  }
  const minus = page.locator("[data-testid='workspace-subtraction-toggle'] button[aria-pressed='true']");
  const difference = page.locator("[data-role='difference-line']");
  while ((await minus.count()) && /^Minus/.test((await minus.first().getAttribute("aria-label")) ?? "") && !(await difference.count()) && Date.now() < until) {
    await page.waitForTimeout(100);
  }
  return {
    background,
    difference_line: await difference.count(),
    minus_pressed: (await minus.count()) > 0 && /^Minus/.test((await minus.first().getAttribute("aria-label")) ?? ""),
    candidates: await page.locator("[data-testid='workspace-candidates']").count(),
    peak_glyphs: await page.locator("[data-role='peak-glyph'] [data-peak-key]").count(),
    custom_index: await page.locator("[data-testid='workspace-custom-index']").count(),
  };
}

export async function readAutosamplerViews(cfg, experimentId, samples) {
  const ui = await reachable(cfg.uiBase, "/");
  if (!ui?.ok) {
    throw new UsageError(`nothing is serving ${cfg.uiBase}`, "--ui reads the workspace in the browser; run `bin/verify seed`");
  }
  await ensureBrowser(cfg);
  const session = await connect(cfg);
  const page = await session.context.newPage();
  const out = new Map();
  try {
    await prepareUser(page, cfg, { username: "verify", onboarding: false });
    const held = new Set(samples.flatMap((s) => s.traceless ?? []));
    for (const s of samples) {
      await settledGoto(page, cfg, WORKSPACE_ROUTE(experimentId, s.id), held);
      out.set(s.id, await readSample(page));
    }
  } finally {
    await page.close().catch(() => undefined);
    await session.browser.close().catch(() => undefined);
  }
  return out;
}
