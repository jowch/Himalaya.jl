import fs from "node:fs";
import path from "node:path";
import { UsageError, intFlag } from "./args.mjs";
import { ensureBrowser, connect } from "./browser.mjs";
import { armSettle } from "./settle.mjs";
import { reachable } from "./http.mjs";

export const STATE_KEY = "himalaya-ui:state";
export const STATE_VERSION = 6;

export function routeUrl(cfg, route) {
  if (/^https?:\/\//.test(route)) return route;
  return new URL(route.startsWith("/") ? route : `/${route}`, cfg.uiBase).toString();
}

export async function leavePage(page) {
  if (page.url() !== "about:blank") await page.goto("about:blank");
}

export async function prepareUser(page, cfg, { username, onboarding }) {
  if (!page.url().startsWith(cfg.uiBase)) await page.goto(cfg.uiBase, { waitUntil: "domcontentloaded" });
  return page.evaluate(
    ({ key, version, username, onboarding }) => {
      const raw = localStorage.getItem(key);
      const parsed = raw ? JSON.parse(raw) : { state: {}, version };
      const state = parsed.state ?? {};
      if (onboarding) {
        parsed.state = { ...state, username: undefined, tutorialSeen: false };
      } else {
        if (state.username && state.tutorialSeen) return "kept";
        parsed.state = { ...state, username: state.username || username, tutorialSeen: true };
      }
      parsed.version = parsed.version ?? version;
      localStorage.setItem(key, JSON.stringify(parsed));
      return onboarding ? "cleared" : "seeded";
    },
    { key: STATE_KEY, version: STATE_VERSION, username, onboarding },
  );
}

export async function withPage(cfg, fn) {
  const session = await connect(cfg);
  try {
    const result = await fn(session.page);
    return result && typeof result === "object" && !Array.isArray(result) ? { ...result, throttle: cfg.throttle } : result;
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
  checkHttpFlags(flags);
  const viewport = parseViewport(flags.viewport);
  const url = routeUrl(cfg, route);
  if (!(await reachable(new URL(url).origin, "/"))) {
    throw new UsageError(`nothing is serving ${new URL(url).origin}`, "run `bin/verify seed`, or `bin/verify doctor` to see what is down");
  }
  const launched = await ensureBrowser(cfg, viewport ?? undefined);
  return withPage(cfg, async (page) => {
    if (viewport) await page.setViewportSize(viewport);
    const user = await prepareUser(page, cfg, { username: flags.user ?? "verify", onboarding: Boolean(flags.onboarding) });
    await leavePage(page);
    const tracker = armSettle(page);
    let response;
    try {
      response = await page.goto(url, { waitUntil: "domcontentloaded" });
    } catch (e) {
      tracker.disarm();
      throw e;
    }
    const settle = await finish(tracker, flags);
    return {
      ...verdict(settle, flags),
      url: page.url(),
      title: await page.title(),
      status: response?.status() ?? null,
      browserStarted: launched.started,
      user,
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
    const nth = intFlag(flags, "nth");
    if (nth >= count) {
      throw new UsageError(`--nth ${flags.nth} is out of range: ${count} element(s) match ${what}`);
    }
    return { loc: loc.nth(nth), what: `${what} nth=${nth}`, count };
  }
  if (count > 1) {
    throw new UsageError(`${count} elements match ${what}`, "add --exact, narrow with --within <css>, or pick one with --nth <i>");
  }
  return { loc, what, count };
}

export function parseExpectHttp(values = []) {
  return values.map((raw) => {
    const m = /^([1-5][0-9x]{2}):(.+)$/.exec(raw);
    if (!m) throw new UsageError(`--expect-http must look like 400:/api/fs/preview, got ${raw}`, "status digits may use x as a wildcard (4xx); the path is a regular expression");
    let path;
    try {
      path = new RegExp(m[2]);
    } catch (e) {
      throw new UsageError(`--expect-http path is not a valid regular expression: ${e.message}`);
    }
    const status = new RegExp(`^${m[1].replace(/x/g, "\\d")}$`);
    return { raw, matches: (e) => status.test(String(e.status)) && path.test(new URL(e.url).pathname) };
  });
}

function checkHttpFlags(flags) {
  if (parseExpectHttp(flags["expect-http"]).length && flags["no-settle"]) {
    throw new UsageError("--expect-http needs the action to settle", "drop --no-settle");
  }
}

const describe = (e) => `${e.method} ${new URL(e.url).pathname} returned ${e.status}`;

export function verdict(settle, flags) {
  const expected = parseExpectHttp(flags["expect-http"]);
  if (!settle) return { ok: true, httpChecked: false };
  for (const e of settle.httpErrors) e.expected = expected.some((x) => x.matches(e));
  const unexpected = settle.httpErrors.filter((e) => !e.expected);
  const unmet = expected.filter((x) => !settle.httpErrors.some((e) => x.matches(e)));
  const problems = [];
  if (!settle.settled) problems.push(`the page did not settle within ${settle.ms} ms`);
  if (unexpected.length) problems.push(`${unexpected.length} request(s) failed, first ${describe(unexpected[0])}`);
  if (settle.settled && unmet.length) problems.push(`no response matched --expect-http ${unmet.map((x) => x.raw).join(", ")}`);
  if (!problems.length) return { ok: true };
  return {
    ok: false,
    error: problems.join("; "),
    hint: unexpected.length
      ? "the page may look fine while the server rejected the change; if this status is the app's normal answer, pass --expect-http <status>:<path regex>"
      : unmet.length && settle.settled ? "the error path you expected did not happen; check the request the action sends" : "raise --timeout, or snapshot to see what is still loading",
  };
}

async function finish(tracker, flags) {
  if (flags["no-settle"]) {
    tracker.disarm();
    return null;
  }
  return tracker.wait({ timeoutMs: intFlag(flags, "timeout", 15000) });
}

async function act(page, flags, action) {
  const tracker = armSettle(page);
  try {
    await action();
  } catch (e) {
    tracker.disarm();
    throw e;
  }
  const settle = await finish(tracker, flags);
  return { ...verdict(settle, flags), url: page.url(), settle };
}

const MODIFIERS = ["Alt", "Control", "ControlOrMeta", "Meta", "Shift"];

export function parseModifiers(value) {
  if (value === undefined) return [];
  const mods = String(value).split(",").map((m) => m.trim()).filter(Boolean);
  const bad = mods.filter((m) => !MODIFIERS.includes(m));
  if (bad.length || mods.length === 0) throw new UsageError(`--modifiers takes a comma list of ${MODIFIERS.join(", ")}, got ${value}`);
  return mods;
}

export async function click(cfg, flags) {
  checkHttpFlags(flags);
  const modifiers = parseModifiers(flags.modifiers);
  return withPage(cfg, async (page) => {
    const target = await resolveOne(page, flags);
    return {
      clicked: target.what,
      ...(modifiers.length ? { modifiers } : {}),
      ...(await act(page, flags, () => target.loc.click({ timeout: 5000, modifiers }))),
    };
  });
}

export async function fill(cfg, flags) {
  checkHttpFlags(flags);
  if (flags.value === undefined) throw new UsageError("fill needs --value <text>");
  return withPage(cfg, async (page) => {
    const target = await resolveOne(page, flags);
    return { filled: target.what, ...(await act(page, flags, () => target.loc.fill(flags.value, { timeout: 5000 }))) };
  });
}

export async function press(cfg, key, flags) {
  checkHttpFlags(flags);
  if (!key) throw new UsageError("press needs a key", "e.g. `bin/verify press Escape` or `bin/verify press Control+K`");
  return withPage(cfg, async (page) => {
    const targeted = flags.role || flags.label || flags.text || flags.testid || flags.selector;
    const target = targeted ? await resolveOne(page, flags) : null;
    const action = target ? () => target.loc.press(key, { timeout: 5000 }) : () => page.keyboard.press(key);
    return { pressed: key, ...(await act(page, flags, action)) };
  });
}

export async function snapshot(cfg, flags) {
  return withPage(cfg, async (page) => {
    const root = flags.selector ? page.locator(flags.selector).first() : page.locator("body");
    const tree = await root.ariaSnapshot({ timeout: 5000 });
    return { url: page.url(), title: await page.title(), snapshot: tree };
  });
}

export async function count(cfg, flags) {
  if (!flags.selector) throw new UsageError("count needs --selector <css>");
  return withPage(cfg, async (page) => {
    const els = page.locator(flags.selector);
    const texts = await els.evaluateAll((nodes) => nodes.slice(0, 20).map((n) => (n.textContent ?? "").trim()));
    return { url: page.url(), selector: flags.selector, count: await els.count(), texts };
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
  checkHttpFlags(flags);
  const timeoutMs = intFlag(flags, "timeout", 15000);
  let pattern = null;
  if (flags.url) {
    try {
      pattern = new RegExp(flags.url);
    } catch (e) {
      throw new UsageError(`--url is not a valid regular expression: ${e.message}`);
    }
  }
  return withPage(cfg, async (page) => {
    const tracker = armSettle(page);
    if (pattern) {
      try {
        await page.waitForURL(pattern, { timeout: timeoutMs, waitUntil: "domcontentloaded" });
      } catch {
        tracker.disarm();
        throw new UsageError(`the page did not reach a URL matching ${flags.url} within ${timeoutMs} ms (now ${page.url()})`, "raise --timeout, or snapshot to see where it is stuck");
      }
    }
    const result = await tracker.wait({ timeoutMs });
    return { ...verdict(result, flags), url: page.url(), ...result };
  });
}
