import { parseArgs, UsageError } from "./args.mjs";
import { loadConfig } from "./config.mjs";
import { doctor } from "./doctor.mjs";
import { seed, stop } from "./seed.mjs";
import { open, click, fill, press, snapshot, screenshot, settle, count } from "./page.mjs";
import { closeBrowser } from "./browser.mjs";
import { api } from "./api.mjs";
import { checkTruth, triage } from "./truth.mjs";
import { deliverLate } from "./late.mjs";

const TARGET = ["role", "name", "label", "text", "testid", "selector", "within", "nth", "timeout"];

export const COMMANDS = {
  doctor: {
    usage: "doctor",
    summary: "Check node, julia, deps, Chromium, demo data, backend, DB, ingest, Vite and its /api proxy.",
    spec: {},
    run: (cfg) => doctor(cfg),
  },
  seed: {
    usage: "seed [--empty] [--reset [--late-trace]] [--dry-run] [--timeout <ms>]",
    summary: "Bring the demo up (generate, backend, ingest, Vite) in the background and wait until doctor passes. No-op when already up. --reset deletes the demo DBs and data first. --late-trace generates the demo with one sample's traces held back in late/, as if the integration pipeline had not written them yet (see deliver-late).",
    spec: { booleans: ["empty", "reset", "dry-run", "late-trace"], strings: ["timeout"] },
    run: (cfg, f) => seed(cfg, f),
  },
  stop: {
    usage: "stop",
    summary: "Stop the demo servers that `seed` started.",
    spec: {},
    run: (cfg) => stop(cfg),
  },
  open: {
    usage: "open <route> [--viewport 1440x900] [--user <name>] [--onboarding] [--no-settle] [--expect-http <status>:<path re>]... [--timeout <ms>]",
    summary: "Open a route in the shared headless browser (started on demand), then wait for it to settle. Seeds a user so onboarding is skipped; --onboarding clears the user so it shows. Fails if the page does not settle. Also fails when any request the page makes while loading returns 4xx/5xx, a missing image included (see click for --expect-http).",
    spec: { booleans: ["onboarding", "no-settle"], strings: ["viewport", "user", "timeout"], lists: ["expect-http"] },
    args: 1,
    run: (cfg, f, [route]) => open(cfg, route, f),
  },
  snapshot: {
    usage: "snapshot [--selector <css>] [--raw]",
    summary: "Accessibility tree (Playwright aria snapshot, YAML) of the page or one element. Use it to find roles and names for click.",
    spec: { booleans: ["raw"], strings: ["selector"] },
    run: (cfg, f) => snapshot(cfg, f),
    raw: (r) => r.snapshot,
  },
  count: {
    usage: "count --selector <css>",
    summary: "Count the elements matching a CSS selector in the open page, right now, with the text of the first 20. Does not wait or settle; poll it to wait for something to appear. Use it for marks the accessibility tree leaves out, such as plot markers.",
    spec: { strings: ["selector"] },
    run: (cfg, f) => count(cfg, f),
  },
  screenshot: {
    usage: "screenshot [<file.png>] [--full] [--selector <css>]",
    summary: "Save a PNG of the viewport, the full page, or one element. Default path scratch/verify/shots/<time>.png.",
    spec: { booleans: ["full"], strings: ["selector"] },
    maxArgs: 1,
    run: (cfg, f, [file]) => screenshot(cfg, file, f),
  },
  click: {
    usage: "click (--role <role> [--name <name>] | --label <l> | --text <t> | --testid <id> | --selector <css>) [--exact] [--within <css>] [--nth <i>] [--modifiers <Control,Meta,Shift,Alt,ControlOrMeta>] [--no-settle] [--expect-http <status>:<path re>]... [--timeout <ms>]",
    summary: "Click exactly one element, holding --modifiers if given (ControlOrMeta is ⌘ on a Mac, Ctrl elsewhere), then wait for the page to settle (requests started by the click included). Fails with a hint when zero or several match, or when the page does not settle within --timeout (default 15000 ms). Also fails when a request the click started returned 4xx/5xx (listed under settle.httpErrors). --expect-http 400:/api/fs/preview (repeatable) allows only matching responses and fails if none arrives. --no-settle skips this check (httpChecked: false).",
    spec: { booleans: ["exact", "no-settle"], strings: [...TARGET, "modifiers"], lists: ["expect-http"] },
    run: (cfg, f) => click(cfg, f),
  },
  fill: {
    usage: "fill <target flags as for click> --value <text> [--no-settle] [--expect-http <status>:<path re>]...",
    summary: "Replace the text of one input, then wait for the page to settle. Fails on 4xx/5xx like click.",
    spec: { booleans: ["exact", "no-settle"], strings: [...TARGET, "value"], lists: ["expect-http"] },
    run: (cfg, f) => fill(cfg, f),
  },
  press: {
    usage: "press <key> [target flags as for click] [--no-settle] [--expect-http <status>:<path re>]...",
    summary: "Press a key (Playwright names: Escape, Enter, ArrowRight, Control+K) on the page or on one element, then settle. Fails on 4xx/5xx like click.",
    spec: { booleans: ["exact", "no-settle"], strings: TARGET, lists: ["expect-http"] },
    args: 1,
    run: (cfg, f, [key]) => press(cfg, key, f),
  },
  "wait-settle": {
    usage: "wait-settle [--url <regex>] [--timeout <ms>] [--expect-http <status>:<path re>]...",
    summary: "Wait until no request is in flight (SSE excluded) and no aria-busy loader or progressbar is visible, for 500 ms. It only sees requests that start after it connects; actions settle on their own. With --url, first wait for a navigation the server triggers later (e.g. after an ingest scan). Fails when it times out, or when a request it saw returned 4xx/5xx (httpErrors at the top level here, not under settle).",
    spec: { strings: ["timeout", "url"], lists: ["expect-http"] },
    run: (cfg, f) => settle(cfg, f),
  },
  close: {
    usage: "close",
    summary: "Close the shared headless browser.",
    spec: {},
    run: (cfg) => closeBrowser(cfg),
  },
  api: {
    usage: "api <METHOD> </api/path> [--data <json>] [--dry-run]",
    summary: "Call the backend and print {ok, status, body}. Writes go through as-is; use --dry-run to see the request first.",
    spec: { booleans: ["dry-run"], strings: ["data"] },
    args: 2,
    run: (cfg, f, [method, p]) => api(cfg, method, p, f),
  },
  "deliver-late": {
    usage: "deliver-late",
    summary: "Move the traces `seed --reset --late-trace` held back into the demo's tot_files, with a fresh mtime, as the integration pipeline would. Then open the experiment so a pass picks them up, and run check-truth.",
    spec: {},
    run: (cfg) => deliverLate(cfg),
  },
  "check-truth": {
    usage: "check-truth [--strict] [--ui [--reload]] [--expect fresh|triaged] [--truth <TRUTH.csv>] [--experiment <id>]",
    summary: "Compare what the app shows with the demo TRUTH.csv: per-frame auto calls, frame counts, and the corpus status (read from GET /api/samples, the data the corpus renders). Fails on any regression against TRUTH's `auto` column; --strict also fails where the auto call misses the generated phases. Also checks an `invariants` block (status matches the confirmed call, representative never a dropped frame, views agree); a broken one fails the exit code. --ui opens each sample in the workspace and reads which frame it opens on and the status it shows; it reuses the shared tab when that is already on the workspace. --reload loads the workspace afresh instead, for a check after a write the open page cannot hear about (a CLI re-index in another process sends no event, so the tab's cache would mix data from before and after it). --expect fresh fails if any frame is dropped; --expect triaged fails if any frame in TRUTH's drop column is still kept. A flared frame's auto call may lose phases but must not add one. A `traces` block fails unless exactly the frames held back in late/ are reported as waiting for a trace (GET /api/experiments/<id>/problems) and no analysis failed; a held frame must have no auto call yet.",
    spec: { booleans: ["strict", "ui", "reload"], strings: ["truth", "experiment", "expect"] },
    run: (cfg, f) => checkTruth(cfg, f),
  },
  triage: {
    usage: "triage [--undo] [--dry-run] [--truth <TRUTH.csv>] [--experiment <id>]",
    summary: "Drop the bad frames listed in TRUTH's `drop` column (PATCH /api/exposures/<id>/status), as a triager would in the workspace. --undo keeps them again, including any of those frames a person dropped before triage ran. Idempotent.",
    spec: { booleans: ["undo", "dry-run"], strings: ["truth", "experiment"] },
    run: (cfg, f) => triage(cfg, f),
  },
};

export function helpText(name) {
  if (name && COMMANDS[name]) {
    const c = COMMANDS[name];
    return `bin/verify ${c.usage}\n\n${c.summary}\n`;
  }
  const lines = Object.entries(COMMANDS).map(([, c]) => `  ${c.usage}\n      ${c.summary}`);
  return [
    "bin/verify <command> [options]   drive and check HimalayaUI on the synthetic demo",
    "",
    "Every command prints one JSON object on stdout and exits 0 on success and 1 on failure; failures carry `error` and a `hint` saying what to do.",
    "Env: DEMO_DIR, DEMO_DB_DIR, DEMO_API_PORT (8080), DEMO_UI_PORT (5173), VERIFY_CDP_PORT (9333), VERIFY_STATE_DIR, VERIFY_CHROMIUM.",
    "",
    "Commands:",
    ...lines,
    "",
    "Typical loop: doctor -> seed -> open /experiments -> snapshot -> click --role link --name ... -> screenshot -> check-truth (after triage, check-truth --ui)",
    "",
  ].join("\n");
}

export async function main(argv, { cfg = loadConfig(), out = process.stdout } = {}) {
  const [name, ...rest] = argv;
  if (!name || name === "help" || name === "--help" || name === "-h") {
    out.write(helpText(rest[0]));
    return 0;
  }
  const cmd = COMMANDS[name];
  if (!cmd) {
    emit(out, { ok: false, error: `unknown command ${name}`, hint: `commands: ${Object.keys(COMMANDS).join(", ")}; run \`bin/verify help\`` });
    return 1;
  }
  if (rest.includes("--help") || rest.includes("-h")) {
    out.write(helpText(name));
    return 0;
  }
  try {
    const { flags, positional } = parseArgs(rest, cmd.spec);
    const want = cmd.args ?? 0;
    const max = cmd.maxArgs ?? want;
    if (positional.length > max) {
      throw new UsageError(`too many arguments: ${positional.slice(max).join(" ")}`, `usage: bin/verify ${cmd.usage}`);
    }
    if (positional.length < want) throw new UsageError("missing argument", `usage: bin/verify ${cmd.usage}`);
    const result = await cmd.run(cfg, flags, positional);
    if (flags.raw && cmd.raw) {
      out.write(cmd.raw(result) + "\n");
      return 0;
    }
    const ok = result.ok !== false;
    emit(out, { ok, command: name, ...result });
    return ok ? 0 : 1;
  } catch (err) {
    if (err instanceof UsageError) {
      emit(out, { ok: false, command: name, error: err.message, ...(err.hint ? { hint: err.hint } : {}) });
      return 1;
    }
    emit(out, { ok: false, command: name, error: err?.message ?? String(err) });
    return 1;
  }
}

function emit(out, obj) {
  out.write(JSON.stringify(obj, null, 2) + "\n");
}
