import { test } from "node:test";
import assert from "node:assert/strict";
import { EventEmitter } from "node:events";
import { armSettle } from "../settle.mjs";
import { findChromium, loadPlaywright } from "../browser.mjs";
import { listen } from "./stub.mjs";

function chromiumOrNull() {
  try {
    return findChromium().path;
  } catch {
    return null;
  }
}
const skip = chromiumOrNull() ? false : "no Chromium available";

class FakePage extends EventEmitter {
  constructor() {
    super();
    this.main = { name: "main" };
  }
  mainFrame() {
    return this.main;
  }
  async evaluate() {
    return { ready: "complete", busy: 0 };
  }
  off(event, fn) {
    return this.removeListener(event, fn);
  }
}

function fakeRequest(url, frame, { navigation = false } = {}) {
  return { url: () => url, resourceType: () => (navigation ? "document" : "fetch"), isNavigationRequest: () => navigation, frame: () => frame };
}

test("a request the old document left open does not hold the next page's settle", async () => {
  const page = new FakePage();
  const tracker = armSettle(page);
  page.emit("request", fakeRequest("http://ui/api/experiments", page.main));
  const nav = fakeRequest("http://ui/experiments/1/workspace/1", page.main, { navigation: true });
  page.emit("request", nav);
  page.emit("requestfinished", nav);
  const settle = await tracker.wait({ quietMs: 50, timeoutMs: 1000 });
  assert.equal(settle.settled, true, `still in flight: ${settle.inflight?.join(" ")}`);
});

test("requests the new page starts after its navigation request still hold the settle", async () => {
  const page = new FakePage();
  const tracker = armSettle(page);
  page.emit("request", fakeRequest("http://ui/api/experiments", page.main));
  const nav = fakeRequest("http://ui/experiments/1/workspace/1", page.main, { navigation: true });
  page.emit("request", nav);
  page.emit("requestfinished", nav);
  page.emit("request", fakeRequest("http://ui/api/experiments/1", page.main));
  const settle = await tracker.wait({ quietMs: 50, timeoutMs: 300 });
  assert.equal(settle.settled, false);
  assert.deepEqual(settle.inflight, ["http://ui/api/experiments/1"]);
});

test("a subframe navigation does not drop the page's requests", async () => {
  const page = new FakePage();
  const tracker = armSettle(page);
  page.emit("request", fakeRequest("http://ui/api/experiments", page.main));
  page.emit("request", fakeRequest("http://ui/embed", { name: "child" }, { navigation: true }));
  const settle = await tracker.wait({ quietMs: 50, timeoutMs: 300 });
  assert.equal(settle.settled, false);
  assert.ok(settle.inflight.includes("http://ui/api/experiments"));
});

test("a request the current page leaves open still holds its settle", { skip, timeout: 30000 }, async () => {
  const site = await listen((req, res) => {
    if (req.url === "/hang") return;
    res.setHeader("Content-Type", "text/html");
    res.end("<!doctype html><title>Stub</title><script>fetch('/hang')</script>");
  });
  const { chromium } = loadPlaywright();
  const browser = await chromium.launch({ executablePath: chromiumOrNull(), args: ["--no-sandbox"] });
  try {
    const page = await browser.newPage();
    const tracker = armSettle(page);
    await page.goto(`http://127.0.0.1:${site.port}/a`, { waitUntil: "domcontentloaded" });
    const settle = await tracker.wait({ timeoutMs: 2000 });
    assert.equal(settle.settled, false);
    assert.deepEqual(settle.inflight, [`http://127.0.0.1:${site.port}/hang`]);
  } finally {
    await browser.close();
    await site.close();
  }
});
