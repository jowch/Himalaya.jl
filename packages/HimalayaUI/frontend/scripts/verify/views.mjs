import { UsageError } from "./args.mjs";
import { ensureBrowser, connect } from "./browser.mjs";
import { armSettle } from "./settle.mjs";
import { reachable } from "./http.mjs";
import { leavePage, prepareUser, routeUrl } from "./page.mjs";

export function statusFromText(text) {
  const t = String(text ?? "").replace(/\s+/g, " ").trim();
  if (!t) return null;
  if (/^(Not indexed|Unreviewed)$|^Index\b/.test(t)) return "unreviewed";
  if (/^Form factor$/i.test(t)) return "form_factor";
  return "phase";
}

export function heldTrace404(e, heldIds) {
  const m = /^\/api\/exposures\/(\d+)\/trace$/.exec(new URL(e.url).pathname);
  return e.method === "GET" && e.status === 404 && m !== null && heldIds.has(Number(m[1]));
}

function failedRequests(settle, where, heldIds = new Set()) {
  const failed = settle.httpErrors.filter((e) => !heldTrace404(e, heldIds));
  if (!failed.length) return;
  const list = failed.map((e) => `${e.method} ${new URL(e.url).pathname} ${e.status}`).join(", ");
  throw new Error(`${where} loaded with failed requests (${list}); the view it shows can't be trusted`);
}

export function frameFromLabels(labels, frameOrder) {
  const numbers = labels.map((l) => /^Image (\d+)(?:,|$)/.exec(l)).filter(Boolean).map((m) => Number(m[1]));
  if (numbers.length !== 1) return null;
  return frameOrder[numbers[0] - 1] ?? null;
}

export const WORKSPACE_ROUTE = (experimentId, id) => `/experiments/${experimentId}/workspace/${id}`;

async function settledGoto(page, cfg, route) {
  await leavePage(page);
  const tracker = armSettle(page);
  try {
    await page.goto(routeUrl(cfg, route), { waitUntil: "domcontentloaded" });
  } catch (e) {
    tracker.disarm();
    throw e;
  }
  const settle = await tracker.wait({ timeoutMs: 20000 });
  if (!settle.settled) throw new Error(`${route} did not settle within 20000 ms (busy ${settle.busy}, in flight: ${settle.inflight.join(" ")})`);
  failedRequests(settle, route);
}

export function workspaceStart(currentUrl, uiBase, experimentId, firstSampleId, { reload = false } = {}) {
  if (reload) return { reuse: false, route: WORKSPACE_ROUTE(experimentId, firstSampleId) };
  let url;
  try {
    url = new URL(currentUrl);
  } catch {
    url = null;
  }
  const prefix = `/experiments/${experimentId}/workspace`;
  if (url && url.origin === new URL(uiBase).origin && (url.pathname === prefix || url.pathname.startsWith(`${prefix}/`))) {
    return { reuse: true, route: null };
  }
  return { reuse: false, route: WORKSPACE_ROUTE(experimentId, firstSampleId) };
}

async function readWorkspaceRows(page) {
  return page.locator("[data-testid='workspace-samples'] [data-sample-id]").evaluateAll((els) =>
    els.map((el) => {
      const mark = el.querySelector("[data-testid='workspace-status']");
      return {
        id: Number(el.getAttribute("data-sample-id")),
        text: mark ? (mark.getAttribute("aria-label") ?? mark.textContent ?? "").trim() : null,
      };
    }),
  );
}

export function keptFromLabels(labels) {
  const frames = labels.filter((l) => /^Image \d+(?:,|$)/.test(l));
  if (!frames.length) return null;
  return `${frames.filter((l) => !/, dropped$/.test(l)).length}/${frames.length}`;
}

async function workspaceKept(page) {
  const labels = await page
    .locator("[data-testid='workspace-exposures'] [aria-label^='Image ']")
    .evaluateAll((els) => els.map((e) => e.getAttribute("aria-label") ?? ""));
  return keptFromLabels(labels);
}

async function workspaceRail(page) {
  const rows = await page.locator("[data-testid='workspace-candidates'] button[aria-pressed]").evaluateAll((els) =>
    els.map((el) => ({
      phase: (el.querySelector("[data-testid='phase-chip']")?.textContent ?? "").trim(),
      on: el.getAttribute("aria-pressed") === "true",
    })),
  );
  const ff = page.locator("[data-testid='workspace-form-factor']");
  const formFactor = (await ff.count()) ? (await ff.getAttribute("aria-pressed")) === "true" : false;
  return { ticked: rows.filter((r) => r.on).map((r) => r.phase), candidates: rows.map((r) => r.phase), form_factor: formFactor };
}

async function currentWorkspaceFrame(page, frameOrder) {
  const strip = "[data-testid='workspace-exposures']";
  const selected = await page.locator(`${strip} [data-selected='true'] button[aria-label^='Image']`)
    .evaluateAll((els) => els.map((e) => e.getAttribute("aria-label") ?? ""));
  if (selected.length === 1) return frameFromLabels(selected, frameOrder);
  const current = page.locator(`${strip} [aria-current='true']`);
  const labels = await current.evaluateAll((els) => els.map((e) => e.getAttribute("aria-label") ?? e.textContent ?? ""));
  return frameFromLabels(labels, frameOrder);
}

export function openSampleBy(currentUrl, experimentId, sampleId) {
  let path = null;
  try {
    path = new URL(currentUrl).pathname;
  } catch {
    path = null;
  }
  return path === WORKSPACE_ROUTE(experimentId, sampleId) ? { reopen: true } : { click: true };
}

async function clickSample(page, id, held = new Set()) {
  const tracker = armSettle(page);
  await page.locator(`[data-testid='workspace-samples'] [data-sample-id='${id}']`).click();
  const settle = await tracker.wait({ timeoutMs: 20000 });
  if (!settle.settled) throw new Error(`the workspace did not settle on sample ${id} within 20000 ms (busy ${settle.busy}, in flight: ${settle.inflight.join(" ")})`);
  failedRequests(settle, `the workspace on sample ${id}`, held);
}

async function textOrNull(page, selector) {
  const el = page.locator(selector);
  return (await el.count()) ? ((await el.first().textContent()) ?? "").trim() : null;
}

async function plotMessage(page, timeoutMs = 5000) {
  const until = Date.now() + timeoutMs;
  for (;;) {
    const text = await textOrNull(page, "[data-testid='workspace-trace-empty']");
    if (text !== "" || Date.now() > until) return text;
    await page.waitForTimeout(100);
  }
}

export async function readWorkspace(page, cfg, experimentId, samples, { reload = false } = {}) {
  const out = new Map();
  if (!samples.length) return out;
  const start = workspaceStart(page.url(), cfg.uiBase, experimentId, samples[0].id, { reload });
  if (!start.reuse) await settledGoto(page, cfg, start.route);
  const filter = page.getByLabel("Filter samples");
  if ((await filter.count()) && (await filter.inputValue())) {
    await filter.fill("");
  }
  const rows = await readWorkspaceRows(page);
  const texts = new Map(rows.map((r) => [r.id, r.text]));
  for (const s of samples) {
    const row = page.locator(`[data-testid='workspace-samples'] [data-sample-id='${s.id}']`);
    let frame = null;
    let kept = null;
    let rail = null;
    let plot = null;
    let notice = null;
    let members = [];
    if (await row.count()) {
      const how = openSampleBy(page.url(), experimentId, s.id);
      if (how.reopen) {
        kept = await workspaceKept(page);
        const away = rows.find((r) => r.id !== s.id);
        if (away) {
          await clickSample(page, away.id, new Set(samples.find((x) => x.id === away.id)?.heldIds ?? []));
          await clickSample(page, s.id, new Set(s.heldIds ?? []));
        } else {
          await settledGoto(page, cfg, WORKSPACE_ROUTE(experimentId, s.id));
        }
      } else {
        await clickSample(page, s.id, new Set(s.heldIds ?? []));
        kept = await workspaceKept(page);
      }
      frame = await currentWorkspaceFrame(page, s.frameOrder);
      members = await page.locator("[data-testid='workspace-members'] [data-member]").evaluateAll((els) => els.map((e) => Number(e.getAttribute("data-member"))));
      if (members.length) frame = "combined";
      rail = await workspaceRail(page);
      plot = await plotMessage(page);
      notice = await textOrNull(page, "[data-testid='workspace-ingest-status']");
    }
    out.set(s.id, { workspace_text: texts.get(s.id) ?? null, workspace: frame, workspace_kept: kept, workspace_rail: rail, workspace_reused: start.reuse, workspace_plot: plot, workspace_notice: notice, workspace_members: members });
  }
  return out;
}

export async function readViews(cfg, experimentId, samples, { reload = false } = {}) {
  if (!(await reachable(cfg.uiBase, "/"))) {
    throw new UsageError(`nothing is serving ${cfg.uiBase}`, "--ui reads the workspace in the browser; run `bin/verify seed`");
  }
  await ensureBrowser(cfg);
  const session = await connect(cfg);
  const out = new Map();
  try {
    await prepareUser(session.page, cfg, { username: "verify", onboarding: false });
    const workspace = await readWorkspace(session.page, cfg, experimentId, samples, { reload });
    for (const s of samples) out.set(s.id, workspace.get(s.id) ?? {});
  } finally {
    await session.browser.close().catch(() => undefined);
  }
  return out;
}
