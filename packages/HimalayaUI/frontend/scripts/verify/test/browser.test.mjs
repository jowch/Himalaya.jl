import { test, before, after } from "node:test";
import assert from "node:assert/strict";
import http from "node:http";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { main } from "../cli.mjs";
import { loadConfig } from "../config.mjs";
import { findChromium, closeBrowser } from "../browser.mjs";

const PAGE = `<!doctype html><html><head><title>Stub</title></head><body>
<main>
  <h1 id="h">Start</h1>
  <button onclick="busy()">Load</button>
  <button>Twin</button><button>Twin</button>
  <button onclick="setTimeout(() => history.pushState({}, '', '/done'), 800)">Later</button>
  <button onclick="fetch('/silent').then((r) => r.text()).then((t) => { document.getElementById('h').textContent = t; })">Silent</button>
  <button onclick="fetch('/hang')">Hang</button>
  <button onclick="fetch('/api/reject', { method: 'PATCH' })">Rejected</button>
  <button onclick="fetch('/api/reject', { method: 'PATCH' }); fetch('/hang')">StuckReject</button>
  <button onclick="setTimeout(() => { fetch('/api/reject', { method: 'PATCH' }); history.pushState({}, '', '/rejected'); }, 800)">DelayedReject</button>
  <label>Probe <input id="q" oninput="fetch('/api/reject?q=' + encodeURIComponent(this.value))"></label>
  <label>Path <input id="p"></label>
  <section id="out" aria-busy="false"></section>
</main>
<script>
  function busy() {
    const out = document.getElementById("out");
    out.setAttribute("aria-busy", "true");
    fetch("/slow").then((r) => r.text()).then((t) => {
      out.setAttribute("aria-busy", "false");
      document.getElementById("h").textContent = t;
    });
  }
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") document.getElementById("h").textContent = "Escaped";
    if (e.key === "Delete") fetch("/api/reject", { method: "DELETE" });
  });
  document.getElementById("p").addEventListener("input", (e) => {
    document.getElementById("h").textContent = "Typed " + e.target.value;
  });
</script></body></html>`;

let server;
let rejects = 0;
let cfg;
let out;
function chromiumOrNull() {
  try {
    return findChromium().path;
  } catch {
    return null;
  }
}
const skip = chromiumOrNull() ? false : "no Chromium available";

function sink() {
  const chunks = [];
  return { write: (s) => chunks.push(s), json: () => JSON.parse(chunks.join("")), text: () => chunks.join("") };
}

async function run(...argv) {
  out = sink();
  const code = await main(argv, { cfg, out });
  return { code, r: argv.includes("--raw") ? out.text() : out.json() };
}

before(async () => {
  server = http.createServer((req, res) => {
    if (req.url === "/slow") {
      setTimeout(() => res.end("Loaded"), 700);
      return;
    }
    if (req.url === "/silent") {
      setTimeout(() => res.end("Loaded silently"), 1500);
      return;
    }
    if (req.url === "/hang") return;
    if (req.url === "/broken") {
      res.setHeader("Content-Type", "text/html");
      res.end('<!doctype html><title>Broken</title><img src="/api/missing.png" alt="">');
      return;
    }
    if (req.url === "/api/missing.png") {
      res.statusCode = 404;
      res.end();
      return;
    }
    if (req.url.startsWith("/api/reject")) {
      rejects++;
      res.statusCode = 400;
      res.end('{"error":"no patchable fields"}');
      return;
    }
    res.setHeader("Content-Type", "text/html");
    res.end(PAGE);
  });
  await new Promise((r) => server.listen(0, "127.0.0.1", r));
  const port = server.address().port;
  const state = fs.mkdtempSync(path.join(os.tmpdir(), "verify-browser-"));
  cfg = loadConfig({ DEMO_UI_PORT: String(port), DEMO_API_PORT: "1", VERIFY_CDP_PORT: String(20000 + (process.pid % 20000)), VERIFY_STATE_DIR: state });
});

after(async () => {
  if (cfg) await closeBrowser(cfg);
  server?.closeAllConnections();
  server?.close();
});

test("open, snapshot, click, fill, press, wait-settle, screenshot, close", { skip, timeout: 60000 }, async () => {
  let { code, r } = await run("open", "/");
  assert.equal(code, 0, JSON.stringify(r));
  assert.equal(r.title, "Stub");
  assert.equal(r.browserStarted, true);
  assert.equal(r.user, "seeded");

  ({ r } = await run("open", "/", "--onboarding"));
  assert.equal(r.user, "cleared");
  ({ r } = await run("open", "/"));
  assert.equal(r.user, "seeded");
  ({ r } = await run("open", "/"));
  assert.equal(r.user, "kept");

  ({ r } = await run("snapshot", "--raw"));
  assert.match(r, /heading "Start" \[level=1\]/);
  assert.match(r, /button "Load"/);

  ({ code, r } = await run("click", "--role", "button", "--name", "Load"));
  assert.equal(code, 0, JSON.stringify(r));
  assert.equal(r.settle.settled, true);
  assert.ok(r.settle.ms >= 600, `settled after ${r.settle.ms} ms, before the slow request finished`);
  ({ r } = await run("snapshot", "--selector", "h1"));
  assert.match(r.snapshot, /Loaded/);

  ({ code, r } = await run("click", "--role", "button", "--name", "Silent"));
  assert.equal(code, 0, JSON.stringify(r));
  assert.ok(r.settle.ms >= 1400, `settled after ${r.settle.ms} ms with a 1500 ms request and no busy marker`);
  ({ r } = await run("snapshot", "--selector", "h1"));
  assert.match(r.snapshot, /Loaded silently/);

  ({ code, r } = await run("click", "--role", "button", "--name", "Hang", "--timeout", "1000"));
  assert.equal(code, 1, "an action that never settles must fail");
  assert.equal(r.settle.settled, false);
  assert.match(r.settle.inflight[0], /\/hang$/);

  ({ code, r } = await run("click", "--role", "button", "--name", "Rejected", "--exact"));
  assert.equal(code, 1, "an action whose request the server rejected must fail");
  assert.equal(r.settle.settled, true);
  assert.deepEqual(r.settle.httpErrors.map((e) => [e.method, new URL(e.url).pathname, e.status]), [["PATCH", "/api/reject", 400]]);
  assert.match(r.error, /PATCH \/api\/reject returned 400/);
  ({ code, r } = await run("click", "--role", "button", "--name", "Rejected", "--exact", "--expect-http", "400:^/api/reject$"));
  assert.equal(code, 0, JSON.stringify(r));
  assert.deepEqual(r.settle.httpErrors.map((e) => [e.status, e.expected]), [[400, true]]);
  ({ code, r } = await run("click", "--role", "button", "--name", "Rejected", "--exact", "--expect-http", "4xx:/api/other"));
  assert.equal(code, 1, "an allowance for another path must not hide this 400");
  assert.match(r.error, /PATCH \/api\/reject returned 400/);
  assert.match(r.error, /no response matched --expect-http 4xx:\/api\/other/);
  ({ code, r } = await run("click", "--role", "button", "--name", "Load", "--expect-http", "400:/api/reject"));
  assert.equal(code, 1, "an expected error that never arrives must fail");
  assert.match(r.error, /no response matched --expect-http 400:\/api\/reject/);
  ({ code, r } = await run("click", "--role", "button", "--name", "Load"));
  assert.equal(code, 0, JSON.stringify(r));
  assert.deepEqual(r.settle.httpErrors, []);
  let before = rejects;
  ({ code, r } = await run("click", "--role", "button", "--name", "Rejected", "--exact", "--expect-http", "nope"));
  assert.equal(code, 1);
  assert.match(r.error, /--expect-http must look like/);
  ({ code, r } = await run("click", "--role", "button", "--name", "Rejected", "--exact", "--expect-http", "400:/api/reject", "--no-settle"));
  assert.equal(code, 1);
  assert.match(r.error, /needs the action to settle/);
  ({ code } = await run("fill", "--label", "Probe", "--value", "z", "--expect-http", "400:("));
  assert.equal(code, 1);
  ({ code } = await run("press", "Delete", "--expect-http", "4:"));
  assert.equal(code, 1);
  await new Promise((r) => setTimeout(r, 300));
  assert.equal(rejects, before, "a malformed --expect-http must be refused before the action runs");

  ({ code, r } = await run("click", "--role", "button", "--name", "StuckReject", "--timeout", "1000"));
  assert.equal(code, 1);
  assert.equal(r.settle.settled, false);
  assert.match(r.error, /did not settle/);
  assert.match(r.error, /PATCH \/api\/reject returned 400/);

  ({ code, r } = await run("click", "--role", "button", "--name", "Rejected", "--exact", "--no-settle"));
  assert.equal(code, 0);
  assert.equal(r.httpChecked, false, "--no-settle must say it did not check responses");

  ({ code, r } = await run("click", "--role", "button", "--name", "DelayedReject", "--no-settle"));
  assert.equal(code, 0);
  ({ code, r } = await run("wait-settle", "--url", "/rejected$", "--timeout", "5000"));
  assert.equal(code, 1, "wait-settle must fail on a 4xx it saw");
  assert.deepEqual(r.httpErrors.map((e) => e.status), [400]);
  ({ code } = await run("open", "/"));
  assert.equal(code, 0);

  ({ code, r } = await run("fill", "--label", "Probe", "--value", "x"));
  assert.equal(code, 1, "fill must fail on a 4xx it triggered");
  assert.match(r.error, /GET \/api\/reject returned 400/);
  ({ code } = await run("fill", "--label", "Probe", "--value", "y", "--expect-http", "400:/api/reject"));
  assert.equal(code, 0);

  ({ code, r } = await run("press", "Delete"));
  assert.equal(code, 1, "press must fail on a 4xx it triggered");
  assert.match(r.error, /DELETE \/api\/reject returned 400/);

  ({ code, r } = await run("open", "/broken"));
  assert.equal(code, 1, "open must fail when the page loads a missing resource");
  assert.match(r.error, /GET \/api\/missing\.png returned 404/);
  ({ code } = await run("open", "/"));
  assert.equal(code, 0);

  ({ code, r } = await run("click", "--role", "button", "--name", "Twin"));
  assert.equal(code, 1);
  assert.match(r.error, /2 elements match/);
  assert.match(r.hint, /--nth/);
  ({ code } = await run("click", "--role", "button", "--name", "Twin", "--nth", "1"));
  assert.equal(code, 0);

  ({ code, r } = await run("click", "--role", "link", "--name", "Nope"));
  assert.equal(code, 1);
  assert.match(r.hint, /snapshot/);

  ({ code } = await run("fill", "--label", "Path", "--value", "abc"));
  assert.equal(code, 0);
  ({ r } = await run("snapshot", "--selector", "h1"));
  assert.match(r.snapshot, /Typed abc/);

  ({ code } = await run("press", "Escape"));
  assert.equal(code, 0);
  ({ r } = await run("snapshot", "--selector", "h1"));
  assert.match(r.snapshot, /Escaped/);

  ({ code, r } = await run("wait-settle", "--timeout", "5000"));
  assert.equal(r.settled, true);

  ({ code } = await run("click", "--role", "button", "--name", "Later", "--no-settle"));
  assert.equal(code, 0);
  ({ code, r } = await run("wait-settle", "--url", "/done$", "--timeout", "5000"));
  assert.equal(code, 0, JSON.stringify(r));
  assert.match(r.url, /\/done$/);
  ({ code, r } = await run("wait-settle", "--url", "(", "--timeout", "100"));
  assert.equal(code, 1);
  assert.match(r.error, /not a valid regular expression/);

  const shot = path.join(cfg.stateDir, "s.png");
  ({ code, r } = await run("screenshot", shot));
  assert.equal(code, 0);
  assert.ok(fs.statSync(shot).size > 1000);

  ({ code, r } = await run("close"));
  assert.equal(code, 0, JSON.stringify(r));
  assert.equal(r.closed, true);
  ({ code, r } = await run("snapshot"));
  assert.equal(code, 1);
});
