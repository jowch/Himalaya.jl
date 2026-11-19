import fs from "node:fs";
import path from "node:path";
import { UsageError } from "./args.mjs";
import { ensureBrowser, connect } from "./browser.mjs";
import { waitSettle } from "./settle.mjs";
import { reachable } from "./http.mjs";

export const STATE_KEY = "himalaya-ui:state";
export const STATE_VERSION = 5;

export function routeUrl(cfg, route) {
  if (/^https?:\/\//.test(route)) return route;
  return new URL(route.startsWith("/") ? route : `/${route}`, cfg.uiBase).toString();
}

async function seedUser(page, cfg, username) {
  if (!page.url().startsWith(cfg.uiBase)) await page.goto(cfg.uiBase, { waitUntil: "domcontentloaded" });
  return page.evaluate(
    ({ key, version, username }) => {
      const raw = localStorage.getItem(key);
      const parsed = raw ? JSON.parse(raw) : { state: {}, version };
      if (parsed.state?.username && parsed.state?.tutorialSeen) return false;
      parsed.state = { ...parsed.state, username: parsed.state?.username || username, tutorialSeen: true };
      parsed.version = parsed.version ?? version;
      localStorage.setItem(key, JSON.stringify(parsed));
      return true;
    },
    { key: STATE_KEY, version: STATE_VERSION, username },
  );
}

export async function withPage(cfg, fn) {
  const session = await connect(cfg);
  try {
    return await fn(session.page);
  } finally {
    await session.browser.close().catch(() => undefined);
  }
}

function parseViewport(value) {
  if (!value) return null;
  const m = /^(\d+)x(\d+)$/.exec(value);
  if (!m) throw new UsageError(`--viewport must look like 1440x900, got ${value}`);
  return { width: Number(m[1]), height: Number(m[2]) };
}

export async function open(cfg, route, flags) {
  const viewport = parseViewport(flags.viewport);
  const url = routeUrl(cfg, route);
  if (!(await reachable(new URL(url).origin, "/"))) {
    throw new UsageError(`nothing is serving ${new URL(url).origin}`, "run `bin/verify seed`, or `bin/verify doctor` to see what is down");
  }
  const launched = await ensureBrowser(cfg, viewport ?? undefined);
  return withPage(cfg, async (page) => {
    if (viewport) await page.setViewportSize(viewport);
    const seeded = flags.onboarding ? false : await seedUser(page, cfg, flags.user ?? "verify");
    const response = await page.goto(url, { waitUntil: "domcontentloaded" });
    const settle = flags["no-settle"] ? null : await waitSettle(page);
    return {
      url: page.url(),
      title: await page.title(),
      status: response?.status() ?? null,
      browserStarted: launched.started,
      seededUser: seeded,
      settle,
    };
  });
}

export function locatorFor(page, flags) {
  const exact = Boolean(flags.exact);
  let loc;
  let what;
  if (flags.role) {
    loc = page.getByRole(flags.role, flags.name !== undefined ? { name: flags.name, exact } : {});
    what = `role=${flags.role}${flags.name !== undefined ? ` name=${JSON.stringify(flags.name)}` : ""}`;
  } else if (flags.label) {
    loc = page.getByLabel(flags.label, { exact });
    what = `label=${JSON.stringify(flags.label)}`;
  } else if (flags.text) {
    loc = page.getByText(flags.text, { exact });
    what = `text=${JSON.stringify(flags.text)}`;
  } else if (flags.testid) {
    loc = page.getByTestId(flags.testid);
    what = `testid=${flags.testid}`;
  } else if (flags.selector) {
    loc = page.locator(flags.selector);
    what = `selector=${flags.selector}`;
  } else {
    throw new UsageError("no target given", "pass --role <role> [--name <name>], --label, --text, --testid or --selector");
  }
  if (flags.within) {
    loc = page.locator(flags.within).locator(loc);
    what = `${flags.within} >> ${what}`;
  }
  return { loc, what };
}

async function resolveOne(page, flags) {
  const { loc, what } = locatorFor(page, flags);
  const count = await loc.count();
  if (count === 0) {
    throw new UsageError(`nothing matches ${what}`, "run `bin/verify snapshot` to see the roles and names on the page");
  }
  if (flags.nth !== undefined) {
    const nth = Number(flags.nth);
    if (!Number.isInteger(nth) || nth < 0 || nth >= count) {
      throw new UsageError(`--nth ${flags.nth} is out of range: ${count} element(s) match ${what}`);
    }
    return { loc: loc.nth(nth), what: `${what} nth=${nth}`, count };
  }
  if (count > 1) {
    throw new UsageError(`${count} elements match ${what}`, "add --exact, narrow with --within <css>, or pick one with --nth <i>");
  }
  return { loc, what, count };
}

async function afterAction(page, flags) {
  return flags["no-settle"] ? null : waitSettle(page);
}

export async function click(cfg, flags) {
  return withPage(cfg, async (page) => {
    const target = await resolveOne(page, flags);
    await target.loc.click({ timeout: 5000 });
    const settle = await afterAction(page, flags);
    return { clicked: target.what, url: page.url(), settle };
  });
}

export async function fill(cfg, flags) {
  if (flags.value === undefined) throw new UsageError("fill needs --value <text>");
  return withPage(cfg, async (page) => {
    const target = await resolveOne(page, flags);
    await target.loc.fill(flags.value, { timeout: 5000 });
    const settle = await afterAction(page, flags);
    return { filled: target.what, url: page.url(), settle };
  });
}

export async function press(cfg, key, flags) {
  if (!key) throw new UsageError("press needs a key", "e.g. `bin/verify press Escape` or `bin/verify press Control+K`");
  return withPage(cfg, async (page) => {
    if (flags.role || flags.label || flags.text || flags.testid || flags.selector) {
      const target = await resolveOne(page, flags);
      await target.loc.press(key, { timeout: 5000 });
    } else {
      await page.keyboard.press(key);
    }
    const settle = await afterAction(page, flags);
    return { pressed: key, url: page.url(), settle };
  });
}

export async function snapshot(cfg, flags) {
  return withPage(cfg, async (page) => {
    const root = flags.selector ? page.locator(flags.selector).first() : page.locator("body");
    const tree = await root.ariaSnapshot({ timeout: 5000 });
    return { url: page.url(), title: await page.title(), snapshot: tree };
  });
}

export async function screenshot(cfg, file, flags) {
  return withPage(cfg, async (page) => {
    const out = path.resolve(file ?? path.join(cfg.stateDir, "shots", `${new Date().toISOString().replace(/[:.]/g, "-")}.png`));
    fs.mkdirSync(path.dirname(out), { recursive: true });
    if (flags.selector) await page.locator(flags.selector).first().screenshot({ path: out });
    else await page.screenshot({ path: out, fullPage: Boolean(flags.full) });
    return { path: out, url: page.url() };
  });
}

export async function settle(cfg, flags) {
  const timeoutMs = flags.timeout ? Number(flags.timeout) : 15000;
  if (!Number.isFinite(timeoutMs) || timeoutMs <= 0) throw new UsageError(`--timeout must be a positive number of ms, got ${flags.timeout}`);
  let pattern = null;
  if (flags.url) {
    try {
      pattern = new RegExp(flags.url);
    } catch (e) {
      throw new UsageError(`--url is not a valid regular expression: ${e.message}`);
    }
  }
  return withPage(cfg, async (page) => {
    if (pattern) {
      try {
        await page.waitForURL(pattern, { timeout: timeoutMs, waitUntil: "domcontentloaded" });
      } catch {
        throw new UsageError(`the page did not reach a URL matching ${flags.url} within ${timeoutMs} ms (now ${page.url()})`, "raise --timeout, or snapshot to see where it is stuck");
      }
    }
    const result = await waitSettle(page, { timeoutMs });
    return { url: page.url(), ...result, ok: result.settled };
  });
}
