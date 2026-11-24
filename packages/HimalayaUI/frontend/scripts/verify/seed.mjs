import fs from "node:fs";
import path from "node:path";
import { spawn, spawnSync } from "node:child_process";
import { commandOf, alive, sleep, waitGone } from "./proc.mjs";
import { REPO } from "./config.mjs";
import { UsageError, intFlag } from "./args.mjs";
import { doctor } from "./doctor.mjs";
import { reachable } from "./http.mjs";
import { STALE_LAYOUT_HINT, staleDemoLayout } from "./demo.mjs";
import { profileDir, resetBrowserProfile } from "./browser.mjs";

const DEMO_SCRIPT = "packages/HimalayaUI/scripts/demo.sh";

function demoEnv(cfg, flags = {}) {
  return {
    ...process.env,
    ...(flags["late-trace"] ? { DEMO_DATA_FLAGS: "--late-trace" } : {}),
    DEMO_DIR: cfg.demoDir,
    DEMO_DB_DIR: cfg.demoDbDir,
    DEMO_API_PORT: String(cfg.apiPort),
    DEMO_UI_PORT: String(cfg.uiPort),
  };
}

function pidPath(cfg) {
  return path.join(cfg.stateDir, "demo.pid");
}


export function readDemoPid(cfg) {
  const file = pidPath(cfg);
  if (!fs.existsSync(file)) return null;
  let record;
  try {
    record = JSON.parse(fs.readFileSync(file, "utf8"));
  } catch {
    return { stale: true };
  }
  const cmd = Number.isInteger(record.pid) ? commandOf(record.pid) : null;
  if (!cmd || !cmd.includes("demo.sh")) return { stale: true, pid: record.pid };
  return { stale: false, pid: record.pid };
}

const NEEDED = ["backend", "database", "frontend", "proxy"];

function ready(report, empty) {
  const need = empty ? NEEDED : [...NEEDED, "demo-ingested"];
  return need.every((n) => report.checks.find((c) => c.name === n)?.ok);
}

async function experimentCount(cfg) {
  const r = await reachable(cfg.apiBase, "/api/experiments");
  return r?.ok && Array.isArray(r.body) ? r.body.length : null;
}

export function seedPlan(cfg, flags) {
  const steps = [];
  if (flags["late-trace"]) steps.push(`generate the demo with --late-trace, holding one sample's traces in ${path.join(cfg.demoDir, "late")}`);
  if (flags.reset) {
    steps.push(`stop the demo servers started by bin/verify (pid file ${pidPath(cfg)})`);
    steps.push(`delete ${cfg.demoDbDir} and ${cfg.demoDbDir}-empty (demo DBs)`);
    steps.push(`delete ${cfg.demoDir} if it holds a TRUTH.csv (make demo-clean)`);
    steps.push(`close the bin/verify browser and delete its profile ${profileDir(cfg)}, so no open tab or localStorage still shows the old DB`);
  }
  steps.push(`run ${DEMO_SCRIPT}${flags.empty ? " --empty" : ""} detached (log ${path.join(cfg.stateDir, "demo.log")})`);
  steps.push(
    flags.empty
      ? "wait until doctor reports backend, database, frontend and proxy, on a DB with no experiments"
      : "wait until doctor reports backend, database, frontend, proxy and demo-ingested",
  );
  return steps;
}

export async function seed(cfg, flags) {
  const empty = Boolean(flags.empty);
  const timeoutMs = intFlag(flags, "timeout", 600000);
  if (flags["late-trace"] && (!flags.reset || empty)) {
    throw new UsageError("--late-trace needs --reset (and no --empty), since it only applies when the demo data is generated", "run `bin/verify seed --reset --late-trace`");
  }
  if (flags["dry-run"]) return { dryRun: true, steps: seedPlan(cfg, flags) };
  if (!flags.reset && !empty && staleDemoLayout(cfg)) {
    throw new UsageError(`${cfg.demoDir} has the layout from before the demo had two setups`, STALE_LAYOUT_HINT);
  }
  if (flags.reset) {
    const stopped = await stop(cfg);
    if (stopped.ok === false) return stopped;
    const browser = await resetBrowserProfile(cfg);
    if (browser.ok === false) return browser;
    spawnSync("make", ["demo-clean", `DEMO_DIR=${cfg.demoDir}`, `DEMO_DB_DIR=${cfg.demoDbDir}`], { cwd: REPO, stdio: "ignore" });
  }
  const before = await doctor(cfg);
  if (ready(before, empty)) {
    if (empty) {
      const n = await experimentCount(cfg);
      if (n !== 0) {
        throw new UsageError(
          `seed --empty needs an empty DB, but the app on :${cfg.apiPort} has ${n} experiment(s)`,
          "run the empty instance on other ports: DEMO_API_PORT=8081 DEMO_UI_PORT=5174 VERIFY_CDP_PORT=9334 VERIFY_STATE_DIR=scratch/verify-empty bin/verify seed --empty",
        );
      }
    }
    return { started: false, alreadyUp: true, doctor: summary(before) };
  }
  if (before.checks.find((c) => c.name === "backend")?.ok) {
    throw new UsageError(
      `a backend already answers on :${cfg.apiPort} but the demo is not ready (failing: ${failing(before)})`,
      readDemoPid(cfg)?.stale === false ? "run `bin/verify stop`, then seed again" : `stop it (lsof -ti:${cfg.apiPort} | xargs -r kill) or use DEMO_API_PORT / DEMO_UI_PORT`,
    );
  }
  fs.mkdirSync(cfg.stateDir, { recursive: true });
  const logPath = path.join(cfg.stateDir, "demo.log");
  const log = fs.openSync(logPath, "w");
  const child = spawn("bash", [DEMO_SCRIPT, ...(empty ? ["--empty"] : [])], {
    cwd: REPO,
    env: demoEnv(cfg, flags),
    detached: true,
    stdio: ["ignore", log, log],
  });
  child.unref();
  fs.writeFileSync(pidPath(cfg), JSON.stringify({ pid: child.pid, script: DEMO_SCRIPT }));
  const started = Date.now();
  for (;;) {
    if (!alive(child.pid)) {
      throw new UsageError("the demo script exited before the app came up", `see ${logPath}`);
    }
    const report = await doctor(cfg);
    if (ready(report, empty)) {
      return { started: true, pid: child.pid, log: logPath, ms: Date.now() - started, doctor: summary(report) };
    }
    if (Date.now() - started > timeoutMs) {
      throw new UsageError(`the demo was not ready after ${timeoutMs} ms`, `see ${logPath}; failing checks: ${failing(report)}`);
    }
    await sleep(2000);
  }
}

function signalGroup(pid, signal) {
  try {
    process.kill(-pid, signal);
    return true;
  } catch {
    try {
      process.kill(pid, signal);
      return true;
    } catch {
      return false;
    }
  }
}

export async function stop(cfg) {
  const record = readDemoPid(cfg);
  if (!record) return { stopped: false, reason: "no demo started by bin/verify is recorded" };
  if (record.stale) {
    fs.rmSync(pidPath(cfg), { force: true });
    return { stopped: false, reason: `the recorded pid ${record.pid ?? "?"} is not a demo.sh process any more; removed the stale pid file` };
  }
  signalGroup(record.pid, "SIGTERM");
  let gone = await waitGone(record.pid, 5000);
  if (!gone) {
    signalGroup(record.pid, "SIGKILL");
    gone = await waitGone(record.pid, 2000);
  }
  if (!gone) {
    return { ok: false, stopped: false, pid: record.pid, error: `demo.sh (pid ${record.pid}) is still running after SIGKILL`, hint: `kill it by hand: kill -9 -${record.pid}` };
  }
  fs.rmSync(pidPath(cfg), { force: true });
  return { stopped: true, pid: record.pid };
}

function summary(report) {
  return Object.fromEntries(report.checks.map((c) => [c.name, c.ok]));
}

function failing(report) {
  return report.checks.filter((c) => !c.ok).map((c) => c.name).join(", ") || "none";
}
