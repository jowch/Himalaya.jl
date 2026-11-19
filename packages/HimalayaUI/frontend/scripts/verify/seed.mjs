import fs from "node:fs";
import path from "node:path";
import { spawn, spawnSync } from "node:child_process";
import { REPO } from "./config.mjs";
import { UsageError } from "./args.mjs";
import { doctor } from "./doctor.mjs";

const DEMO_SCRIPT = "packages/HimalayaUI/scripts/demo.sh";

function demoEnv(cfg) {
  return {
    ...process.env,
    DEMO_DIR: cfg.demoDir,
    DEMO_DB_DIR: cfg.demoDbDir,
    DEMO_API_PORT: String(cfg.apiPort),
    DEMO_UI_PORT: String(cfg.uiPort),
  };
}

function pidPath(cfg) {
  return path.join(cfg.stateDir, "demo.pid");
}

function readPid(cfg) {
  const file = pidPath(cfg);
  if (!fs.existsSync(file)) return null;
  const pid = Number(fs.readFileSync(file, "utf8"));
  return alive(pid) ? pid : null;
}

function alive(pid) {
  try {
    process.kill(pid, 0);
    return true;
  } catch {
    return false;
  }
}

function ready(report, empty) {
  const need = empty ? ["backend", "database", "frontend", "proxy"] : ["backend", "database", "demo-ingested", "frontend", "proxy"];
  return need.every((n) => report.checks.find((c) => c.name === n)?.ok);
}

export function seedPlan(cfg, flags) {
  const steps = [];
  if (flags.reset) {
    steps.push(`stop the demo servers started by bin/verify (pid file ${pidPath(cfg)})`);
    steps.push(`delete ${cfg.demoDbDir} and ${cfg.demoDbDir}-empty (demo DBs)`);
    steps.push(`delete ${cfg.demoDir} if it holds a TRUTH.csv (make demo-clean)`);
  }
  steps.push(`run ${DEMO_SCRIPT}${flags.empty ? " --empty" : ""} detached (log ${path.join(cfg.stateDir, "demo.log")})`);
  steps.push("wait until doctor reports backend, database, frontend, proxy" + (flags.empty ? "" : " and demo-ingested"));
  return steps;
}

export async function seed(cfg, flags) {
  const empty = Boolean(flags.empty);
  if (flags["dry-run"]) return { dryRun: true, steps: seedPlan(cfg, flags) };
  if (flags.reset) {
    await stop(cfg);
    spawnSync("make", ["demo-clean", `DEMO_DIR=${cfg.demoDir}`, `DEMO_DB_DIR=${cfg.demoDbDir}`], { cwd: REPO, stdio: "ignore" });
  }
  const before = await doctor(cfg);
  if (!flags.reset && ready(before, empty)) {
    return { started: false, alreadyUp: true, doctor: summary(before) };
  }
  const backend = before.checks.find((c) => c.name === "backend");
  if (backend?.ok) {
    throw new UsageError(
      `a backend already answers on :${cfg.apiPort} but the demo is not ready`,
      readPid(cfg) ? "run `bin/verify stop`, then seed again" : `stop it (lsof -ti:${cfg.apiPort} | xargs -r kill) or use DEMO_API_PORT / DEMO_UI_PORT`,
    );
  }
  fs.mkdirSync(cfg.stateDir, { recursive: true });
  const logPath = path.join(cfg.stateDir, "demo.log");
  const log = fs.openSync(logPath, "w");
  const child = spawn("bash", [DEMO_SCRIPT, ...(empty ? ["--empty"] : [])], {
    cwd: REPO,
    env: demoEnv(cfg),
    detached: true,
    stdio: ["ignore", log, log],
  });
  child.unref();
  fs.writeFileSync(pidPath(cfg), String(child.pid));
  const timeoutMs = Number(flags.timeout ?? 600000);
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
    await new Promise((r) => setTimeout(r, 2000));
  }
}

export async function stop(cfg) {
  const pid = readPid(cfg);
  if (fs.existsSync(pidPath(cfg))) fs.rmSync(pidPath(cfg));
  if (!pid) return { stopped: false, reason: "no demo started by bin/verify is running" };
  try {
    process.kill(-pid, "SIGTERM");
  } catch {
    process.kill(pid, "SIGTERM");
  }
  for (let i = 0; i < 50 && alive(pid); i++) await new Promise((r) => setTimeout(r, 100));
  return { stopped: true, pid };
}

function summary(report) {
  return Object.fromEntries(report.checks.map((c) => [c.name, c.ok]));
}

function failing(report) {
  return report.checks.filter((c) => !c.ok).map((c) => c.name).join(", ") || "none";
}
