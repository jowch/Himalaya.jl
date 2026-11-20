export const BUSY_SELECTOR = "[data-boneyard-overlay], [aria-busy='true'], [role='progressbar']";

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

export function tracked(request) {
  return !["eventsource", "font"].includes(request.resourceType()) && !/^(blob|data):/.test(request.url());
}

export function armSettle(page) {
  const inflight = new Set();
  const onRequest = (r) => {
    if (tracked(r)) inflight.add(r);
  };
  const onDone = (r) => inflight.delete(r);
  page.on("request", onRequest);
  page.on("requestfinished", onDone);
  page.on("requestfailed", onDone);
  const disarm = () => {
    page.off("request", onRequest);
    page.off("requestfinished", onDone);
    page.off("requestfailed", onDone);
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
          return { settled: true, ms: now - started };
        }
        if (now - started > timeoutMs) {
          return { settled: false, ms: now - started, busy: probe.busy, inflight: [...inflight].map((r) => r.url()).slice(0, 5) };
        }
        await sleep(100);
      }
    } finally {
      disarm();
    }
  }

  return { wait, disarm };
}
