export const SKELETON_SELECTOR = "[data-boneyard-overlay], [aria-busy='true'], [role='progressbar']";

export async function waitSettle(page, { quietMs = 500, timeoutMs = 15000 } = {}) {
  const started = Date.now();
  const inflight = new Set();
  const onRequest = (r) => {
    if (r.resourceType() !== "eventsource") inflight.add(r);
  };
  const onDone = (r) => inflight.delete(r);
  page.on("request", onRequest);
  page.on("requestfinished", onDone);
  page.on("requestfailed", onDone);
  try {
    let quietSince = Date.now();
    let lastEntries = -1;
    for (;;) {
      const probe = await page.evaluate((sel) => {
        const entries = performance.getEntriesByType("resource").length;
        const skeletons = [...document.querySelectorAll(sel)].filter((el) => {
          const s = getComputedStyle(el);
          return s.display !== "none" && s.visibility !== "hidden" && el.getClientRects().length > 0;
        }).length;
        return { ready: document.readyState, entries, skeletons };
      }, SKELETON_SELECTOR);
      const busy = probe.ready !== "complete" || probe.skeletons > 0 || inflight.size > 0 || probe.entries !== lastEntries;
      lastEntries = probe.entries;
      const now = Date.now();
      if (busy) quietSince = now;
      if (now - quietSince >= quietMs) {
        await page.evaluate(() => new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(() => r()))));
        return { settled: true, ms: now - started, skeletons: 0 };
      }
      if (now - started > timeoutMs) {
        return { settled: false, ms: now - started, skeletons: probe.skeletons, inflight: inflight.size };
      }
      await new Promise((r) => setTimeout(r, 100));
    }
  } finally {
    page.off("request", onRequest);
    page.off("requestfinished", onDone);
    page.off("requestfailed", onDone);
  }
}
