import { UsageError } from "./args.mjs";
import { ensureBrowser, connect } from "./browser.mjs";
import { armSettle } from "./settle.mjs";
import { reachable } from "./http.mjs";
import { prepareUser, routeUrl } from "./page.mjs";

export const VIEW_ROUTES = {
  loupe: (id) => `/sample/${id}/loupe`,
  focus: (id) => `/sample/${id}`,
};

export function statusFromText(text) {
  const t = String(text ?? "").replace(/\s+/g, " ").trim();
  if (!t) return null;
  if (/^(Not indexed|Unreviewed)$|^Index\b/.test(t)) return "unreviewed";
  if (/^Form factor$/i.test(t)) return "form_factor";
  if (/^No exposures$/.test(t)) return "no_exposures";
  return "phase";
}

async function corpusStatuses(page, cfg, experimentId) {
  const tracker = armSettle(page);
  await page.goto(routeUrl(cfg, `/experiments/${experimentId}/corpus`), { waitUntil: "domcontentloaded" });
  const settle = await tracker.wait({ timeoutMs: 20000 });
  if (!settle.settled) throw new Error(`the corpus did not settle within 20000 ms (busy ${settle.busy}, in flight: ${settle.inflight.join(" ")})`);
  return page.locator("[role='row']").evaluateAll((rows) =>
    rows
      .map((row) => ({
        names: [...row.querySelectorAll("button")].map((b) => (b.textContent ?? "").trim()),
        status: row.querySelector("[data-testid='status-cell']")?.textContent?.trim() ?? null,
        kept: (() => {
          const cell = row.querySelector("[data-testid='kept-cell']");
          const count = cell?.querySelector("[data-role='kept-count']")?.textContent?.trim();
          const total = cell?.querySelector("[data-role='kept-total']")?.textContent?.replace(/[^0-9]/g, "");
          return count && total ? `${count}/${total}` : null;
        })(),
      }))
      .filter((r) => r.status !== null),
  );
}

export function frameFromLabels(labels, frameOrder) {
  const numbers = labels.map((l) => /^Frame (\d+)\b/.exec(l)).filter(Boolean).map((m) => Number(m[1]));
  if (numbers.length !== 1) return null;
  return frameOrder[numbers[0] - 1] ?? null;
}

async function openedFrame(page, cfg, route, frameOrder) {
  await settledGoto(page, cfg, route);
  const pressed = page.getByRole("button", { name: /^Frame \d+/, pressed: true });
  const labels = await pressed.evaluateAll((els) => els.map((e) => e.getAttribute("aria-label") ?? e.textContent ?? ""));
  return frameFromLabels(labels, frameOrder);
}

export const WORKSPACE_ROUTE = (experimentId, id) => `/experiments/${experimentId}/workspace/${id}`;

async function settledGoto(page, cfg, route) {
  const tracker = armSettle(page);
  try {
    await page.goto(routeUrl(cfg, route), { waitUntil: "domcontentloaded" });
  } catch (e) {
    tracker.disarm();
    throw e;
  }
  const settle = await tracker.wait({ timeoutMs: 20000 });
  if (!settle.settled) throw new Error(`${route} did not settle within 20000 ms (busy ${settle.busy}, in flight: ${settle.inflight.join(" ")})`);
}

export function workspaceStart(currentUrl, uiBase, experimentId, firstSampleId) {
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

async function currentWorkspaceFrame(page, frameOrder) {
  const current = page.locator("[data-testid='workspace-exposures'] [aria-current='true']");
  const labels = await current.evaluateAll((els) => els.map((e) => e.getAttribute("aria-label") ?? e.textContent ?? ""));
  return frameFromLabels(labels, frameOrder);
}

export async function readWorkspace(page, cfg, experimentId, samples) {
  const out = new Map();
  if (!samples.length) return out;
  const start = workspaceStart(page.url(), cfg.uiBase, experimentId, samples[0].id);
  if (!start.reuse) await settledGoto(page, cfg, start.route);
  const filter = page.getByLabel("Filter samples");
  if ((await filter.count()) && (await filter.inputValue())) {
    await filter.fill("");
  }
  const texts = new Map((await readWorkspaceRows(page)).map((r) => [r.id, r.text]));
  for (const s of samples) {
    const row = page.locator(`[data-testid='workspace-samples'] [data-sample-id='${s.id}']`);
    let frame = null;
    if (await row.count()) {
      const tracker = armSettle(page);
      await row.click();
      const settle = await tracker.wait({ timeoutMs: 20000 });
      if (!settle.settled) throw new Error(`the workspace did not settle on sample ${s.id} within 20000 ms (busy ${settle.busy}, in flight: ${settle.inflight.join(" ")})`);
      frame = await currentWorkspaceFrame(page, s.frameOrder);
    }
    out.set(s.id, { workspace_text: texts.get(s.id) ?? null, workspace: frame, workspace_reused: start.reuse });
  }
  return out;
}

export async function readViews(cfg, experimentId, samples) {
  if (!(await reachable(cfg.uiBase, "/"))) {
    throw new UsageError(`nothing is serving ${cfg.uiBase}`, "--ui reads the corpus, Loupe, Focus and the workspace in the browser; run `bin/verify seed`");
  }
  await ensureBrowser(cfg);
  const session = await connect(cfg);
  const out = new Map();
  try {
    await prepareUser(session.page, cfg, { username: "verify", onboarding: false });
    const workspace = await readWorkspace(session.page, cfg, experimentId, samples);
    const rows = await corpusStatuses(session.page, cfg, experimentId);
    for (const s of samples) {
      const row = rows.filter((r) => r.names.includes(s.name));
      const views = { corpus_text: row.length === 1 ? row[0].status : null, corpus_kept: row.length === 1 ? row[0].kept : null };
      for (const [view, route] of Object.entries(VIEW_ROUTES)) {
        views[view] = await openedFrame(session.page, cfg, route(s.id), s.frameOrder);
      }
      Object.assign(views, workspace.get(s.id) ?? {});
      out.set(s.id, views);
    }
  } finally {
    await session.browser.close().catch(() => undefined);
  }
  return out;
}
