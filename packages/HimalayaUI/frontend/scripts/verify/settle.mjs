export const BUSY_SELECTOR = "[data-boneyard-overlay], [aria-busy='true'], [role='progressbar']";

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

export function tracked(request) {
  const url = request.url();
  return !["eventsource", "font"].includes(request.resourceType()) && !/^(blob|data):/.test(url) && !/\/favicon\.[a-z]+(\?|$)/.test(url);
}

export function armSettle(page) {
  const inflight = new Set();
  const httpErrors = [];
  let navigating = null;
  const mainNavigation = (r) => r.isNavigationRequest() && r.frame() === page.mainFrame();
  const onRequest = (r) => {
    if (mainNavigation(r)) navigating = r;
    if (tracked(r)) inflight.add(r);
  };
  const onDone = (r) => inflight.delete(r);
  const onFailed = (r) => {
    if (r === navigating) navigating = null;
    inflight.delete(r);
  };
  const onNavigated = (frame) => {
    if (frame !== page.mainFrame() || navigating === null) return;
    navigating = null;
    for (const r of inflight) if (!r.isNavigationRequest()) inflight.delete(r);
  };
  const onResponse = (res) => {
    const r = res.request();
    if (res.status() >= 400 && tracked(r)) httpErrors.push({ method: r.method(), url: r.url(), status: res.status() });
  };
  page.on("request", onRequest);
  page.on("response", onResponse);
  page.on("requestfinished", onDone);
  page.on("requestfailed", onFailed);
  page.on("framenavigated", onNavigated);
  const disarm = () => {
    page.off("request", onRequest);
    page.off("response", onResponse);
    page.off("requestfinished", onDone);
    page.off("requestfailed", onFailed);
    page.off("framenavigated", onNavigated);
  };

  async function wait({ quietMs = 500, timeoutMs = 15000 } = {}) {
    const started = Date.now();
    try {
      let quietSince = Date.now();
      for (;;) {
        const probe = await page.evaluate((sel) => {
          const busy = [...document.querySelectorAll(sel)].filter((el) => {
            const s = getComputedStyle(el);
            return s.display !== "none" && s.visibility !== "hidden" && el.getClientRects().length > 0;
          }).length;
          return { ready: document.readyState, busy };
        }, BUSY_SELECTOR);
        const now = Date.now();
        if (probe.ready !== "complete" || probe.busy > 0 || inflight.size > 0) quietSince = now;
        if (now - quietSince >= quietMs) {
          await page.evaluate(() => new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(() => r()))));
          return { settled: true, ms: now - started, httpErrors };
        }
        if (now - started > timeoutMs) {
          return { settled: false, ms: now - started, busy: probe.busy, inflight: [...inflight].map((r) => r.url()).slice(0, 5), httpErrors };
        }
        await sleep(100);
      }
    } finally {
      disarm();
    }
  }

  return { wait, disarm };
}
