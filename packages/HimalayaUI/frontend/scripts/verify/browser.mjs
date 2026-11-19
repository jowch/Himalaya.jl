import { spawn } from "node:child_process";
import fs from "node:fs";
import path from "node:path";
import { createRequire } from "node:module";
import { FRONTEND } from "./config.mjs";
import { UsageError } from "./args.mjs";

const require = createRequire(path.join(FRONTEND, "package.json"));

export function loadPlaywright() {
  try {
    return require("playwright");
  } catch {
    throw new UsageError(
      "playwright is not installed in the frontend",
      "run `make setup` (or `npm ci` in packages/HimalayaUI/frontend)",
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

export async function browserStatus(cfg) {
  const version = await cdpVersion(cfg);
  return { running: version !== null, cdp: cfg.cdpBase, browser: version?.Browser ?? null };
}

export async function ensureBrowser(cfg, { width = 1440, height = 900 } = {}) {
  if (await cdpVersion(cfg)) return { started: false };
  const exe = findChromium();
  if (!exe.path) {
    throw new UsageError(
      "no Chromium found for Playwright",
      "run `npx playwright install chromium` in packages/HimalayaUI/frontend, or set VERIFY_CHROMIUM=/path/to/chrome",
    );
  }
  fs.mkdirSync(cfg.stateDir, { recursive: true });
  const profile = path.join(cfg.stateDir, "profile");
  const log = fs.openSync(path.join(cfg.stateDir, "browser.log"), "a");
  const child = spawn(
    exe.path,
    [
      "--headless=new",
      "--no-sandbox",
      "--disable-gpu",
      "--no-first-run",
      "--no-default-browser-check",
      `--remote-debugging-port=${cfg.cdpPort}`,
      `--user-data-dir=${profile}`,
      `--window-size=${width},${height}`,
      "about:blank",
    ],
    { detached: true, stdio: ["ignore", log, log] },
  );
  child.unref();
  fs.writeFileSync(pidFile(cfg), String(child.pid));
  for (let i = 0; i < 100; i++) {
    if (await cdpVersion(cfg)) return { started: true, pid: child.pid, chromium: exe };
    await new Promise((r) => setTimeout(r, 100));
  }
  throw new UsageError("Chromium did not open its debugging port", `see ${path.join(cfg.stateDir, "browser.log")}`);
}

export async function connect(cfg) {
  if (!(await cdpVersion(cfg))) {
    throw new UsageError("no browser session is open", "run `bin/verify open <route>` first");
  }
  const { chromium } = loadPlaywright();
  const browser = await chromium.connectOverCDP(cfg.cdpBase);
  const context = browser.contexts()[0] ?? (await browser.newContext());
  const pages = context.pages().filter((p) => !p.url().startsWith("devtools://"));
  const page = pages[pages.length - 1] ?? (await context.newPage());
  return { browser, context, page };
}

export async function closeBrowser(cfg) {
  const file = pidFile(cfg);
  let pid = null;
  if (fs.existsSync(file)) {
    pid = Number(fs.readFileSync(file, "utf8"));
    fs.rmSync(file);
  }
  const wasRunning = (await cdpVersion(cfg)) !== null;
  if (wasRunning) {
    const { chromium } = loadPlaywright();
    const browser = await chromium.connectOverCDP(cfg.cdpBase);
    try {
      const session = await browser.newBrowserCDPSession();
      await session.send("Browser.close");
    } catch {
      if (pid) killQuietly(pid);
    }
  }
  return { closed: wasRunning };
}

function killQuietly(pid) {
  try {
    process.kill(pid);
    return true;
  } catch {
    return false;
  }
}
