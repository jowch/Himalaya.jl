import { spawn } from "node:child_process";
import { applyThrottle } from "./throttle.mjs";
import fs from "node:fs";
import path from "node:path";
import { createRequire } from "node:module";
import { FRONTEND } from "./config.mjs";
import { UsageError } from "./args.mjs";
import { commandOf, sleep, waitGone } from "./proc.mjs";

const require = createRequire(path.join(FRONTEND, "package.json"));

export function loadPlaywright() {
  try {
    return require("playwright");
  } catch {
    throw new UsageError(
      "playwright is not installed in the frontend",
      "run `make dev-setup` (or `npm ci` in packages/HimalayaUI/frontend)",
    );
  }
}

function newestChrome(root) {
  if (!fs.existsSync(root)) return null;
  const dirs = fs
    .readdirSync(root)
    .filter((d) => /^chromium-\d+$/.test(d))
    .sort((a, b) => Number(b.split("-")[1]) - Number(a.split("-")[1]));
  for (const d of dirs) {
    const exe = path.join(root, d, "chrome-linux64", "chrome");
    if (fs.existsSync(exe)) return exe;
  }
  return null;
}

export function findChromium(env = process.env) {
  if (env.VERIFY_CHROMIUM) {
    if (!fs.existsSync(env.VERIFY_CHROMIUM)) return { path: null, source: "VERIFY_CHROMIUM (missing)" };
    return { path: env.VERIFY_CHROMIUM, source: "VERIFY_CHROMIUM" };
  }
  const { chromium } = loadPlaywright();
  const bundled = chromium.executablePath();
  if (fs.existsSync(bundled)) return { path: bundled, source: "playwright" };
  const preinstalled = newestChrome(path.dirname(path.dirname(path.dirname(bundled))));
  if (preinstalled) return { path: preinstalled, source: "preinstalled" };
  return { path: null, source: "none" };
}

async function cdpVersion(cfg) {
  try {
    const res = await fetch(`${cfg.cdpBase}/json/version`, { signal: AbortSignal.timeout(1000) });
    return res.ok ? await res.json() : null;
  } catch {
    return null;
  }
}

function pidFile(cfg) {
  return path.join(cfg.stateDir, "browser.pid");
}

export function profileDir(cfg) {
  return path.join(cfg.stateDir, "profile");
}

function ownPid(cfg) {
  const file = pidFile(cfg);
  if (!fs.existsSync(file)) return null;
  const pid = Number(fs.readFileSync(file, "utf8"));
  const cmd = Number.isInteger(pid) ? commandOf(pid) : null;
  return cmd && cmd.includes(`--user-data-dir=${profileDir(cfg)}`) ? pid : null;
}

export async function browserStatus(cfg) {
  const version = await cdpVersion(cfg);
  return { running: version !== null, ours: version !== null && ownPid(cfg) !== null, cdp: cfg.cdpBase, browser: version?.Browser ?? null };
}

function foreign(cfg) {
  return new UsageError(
    `something not started by bin/verify answers on ${cfg.cdpBase}`,
    "set VERIFY_CDP_PORT to a free port, or stop whatever is listening there",
  );
}

export async function ensureBrowser(cfg, { width = 1440, height = 900 } = {}) {
  if (await cdpVersion(cfg)) {
    if (ownPid(cfg) === null) throw foreign(cfg);
    return { started: false };
  }
  const exe = findChromium();
  if (!exe.path) {
    throw new UsageError(
      "no Chromium found for Playwright",
      "run `npx playwright install chromium` in packages/HimalayaUI/frontend, or set VERIFY_CHROMIUM=/path/to/chrome",
    );
  }
  fs.mkdirSync(cfg.stateDir, { recursive: true });
  const log = fs.openSync(path.join(cfg.stateDir, "browser.log"), "a");
  const child = spawn(
    exe.path,
    [
      "--headless=new",
      "--no-sandbox",
      "--disable-gpu",
      "--no-first-run",
      "--no-default-browser-check",
      "--disable-features=BackForwardCache",
      `--remote-debugging-port=${cfg.cdpPort}`,
      `--user-data-dir=${profileDir(cfg)}`,
      `--window-size=${width},${height}`,
      "about:blank",
    ],
    { detached: true, stdio: ["ignore", log, log] },
  );
  child.unref();
  fs.writeFileSync(pidFile(cfg), String(child.pid));
  for (let i = 0; i < 100; i++) {
    if (await cdpVersion(cfg)) return { started: true, pid: child.pid, chromium: exe };
    await sleep(100);
  }
  throw new UsageError("Chromium did not open its debugging port", `see ${path.join(cfg.stateDir, "browser.log")}`);
}

export async function connect(cfg) {
  if (!(await cdpVersion(cfg))) {
    throw new UsageError("no browser session is open", "run `bin/verify open <route>` first");
  }
  if (ownPid(cfg) === null) throw foreign(cfg);
  const { chromium } = loadPlaywright();
  const browser = await chromium.connectOverCDP(cfg.cdpBase);
  const context = browser.contexts()[0] ?? (await browser.newContext());
  const pages = context.pages().filter((p) => !p.url().startsWith("devtools://"));
  const page = pages[pages.length - 1] ?? (await context.newPage());
  await applyThrottle(context, page, cfg.throttle);
  return { browser, context, page };
}

async function politeClose(cfg) {
  const { chromium } = loadPlaywright();
  const browser = await chromium.connectOverCDP(cfg.cdpBase, { timeout: 3000 });
  const session = await browser.newBrowserCDPSession();
  await session.send("Browser.close").catch(() => undefined);
}

export async function closeBrowser(cfg) {
  const answering = (await cdpVersion(cfg)) !== null;
  const pid = ownPid(cfg);
  if (answering && pid === null) throw foreign(cfg);
  if (pid === null) {
    fs.rmSync(pidFile(cfg), { force: true });
    return { closed: false };
  }
  await Promise.race([politeClose(cfg).catch(() => undefined), sleep(3000)]);
  if (!(await waitGone(pid, 3000))) killQuietly(pid);
  if (!(await waitGone(pid, 2000))) {
    return { ok: false, closed: false, error: `Chromium (pid ${pid}) is still running`, hint: `kill it by hand: kill -9 ${pid}` };
  }
  fs.rmSync(pidFile(cfg), { force: true });
  return { closed: true };
}

export async function resetBrowserProfile(cfg) {
  const closed = await closeBrowser(cfg);
  if (closed.ok === false) return closed;
  fs.rmSync(profileDir(cfg), { recursive: true, force: true, maxRetries: 3 });
  return { ...closed, profile: "deleted" };
}

function killQuietly(pid) {
  try {
    process.kill(pid, "SIGKILL");
    return true;
  } catch {
    return false;
  }
}
