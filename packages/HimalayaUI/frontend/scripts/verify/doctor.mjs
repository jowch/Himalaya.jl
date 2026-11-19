import fs from "node:fs";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { FRONTEND } from "./config.mjs";
import { findChromium, browserStatus } from "./browser.mjs";
import { reachable } from "./http.mjs";
import { readTruth, findDemoExperiment } from "./demo.mjs";

const SEED = "run `bin/verify seed` (starts backend + Vite on the demo data and waits for ingest)";

function check(name, ok, detail, fix) {
  return ok ? { name, ok, detail } : { name, ok, detail, fix };
}

export async function doctor(cfg) {
  const checks = [];
  const nodeMajor = Number(process.versions.node.split(".")[0]);
  checks.push(check("node", nodeMajor >= 20, `v${process.versions.node}`, "install Node >= 20"));

  const jl = spawnSync("julia", ["--version"], { encoding: "utf8" });
  checks.push(
    check("julia", jl.status === 0, jl.status === 0 ? jl.stdout.trim() : "not on PATH", "install Julia (see docs/demo.md, Cloud sessions)"),
  );

  const deps = fs.existsSync(path.join(FRONTEND, "node_modules", "playwright"));
  checks.push(check("frontend-deps", deps, deps ? "node_modules present" : "node_modules missing", "run `make setup`"));

  if (deps) {
    const exe = findChromium();
    checks.push(
      check("chromium", exe.path !== null, exe.path ? `${exe.source}: ${exe.path}` : "no browser binary", "run `npx playwright install chromium` in packages/HimalayaUI/frontend or set VERIFY_CHROMIUM"),
    );
  }

  const truth = readTruth(cfg);
  const images = truth?.every((r) => r.images === "true");
  checks.push(
    check(
      "demo-data",
      Boolean(truth && images),
      truth ? `${truth.length} samples in ${cfg.demoDir}${images ? "" : " (generated with --no-images)"}` : `no TRUTH.csv in ${cfg.demoDir}`,
      "run `make demo-data`",
    ),
  );

  const health = await reachable(cfg.apiBase, "/api/health");
  checks.push(check("backend", Boolean(health?.ok), health ? `${cfg.apiBase} answered ${health.status}` : `${cfg.apiBase} not answering`, SEED));

  if (health?.ok) {
    const exps = await reachable(cfg.apiBase, "/api/experiments");
    const dbOk = Boolean(exps?.ok && Array.isArray(exps.body));
    checks.push(check("database", dbOk, dbOk ? `${exps.body.length} experiment(s)` : `GET /api/experiments returned ${exps?.status}`, "check the backend log in scratch/demo/backend.log"));
    if (dbOk) {
      const exp = await findDemoExperiment(cfg);
      const done = exp?.ingest_status === "complete";
      checks.push(
        check("demo-ingested", done, exp ? `experiment ${exp.id} is ${exp.ingest_status}` : "demo experiment not in the DB", exp ? "wait for the scan, or `bin/verify seed`" : SEED),
      );
    }
  }

  const ui = await reachable(cfg.uiBase, "/");
  const served = Boolean(ui?.ok && typeof ui.body === "string" && ui.body.includes('id="app"'));
  checks.push(check("frontend", served, ui ? `${cfg.uiBase} answered ${ui.status}` : `${cfg.uiBase} not answering`, SEED));
  if (served) {
    const proxied = await reachable(cfg.uiBase, "/api/health");
    checks.push(
      check("proxy", Boolean(proxied?.ok), proxied ? `${cfg.uiBase}/api answered ${proxied.status}` : "no answer", `restart Vite with VITE_API_PORT=${cfg.apiPort}`),
    );
  }

  const browser = await browserStatus(cfg);
  const ok = checks.every((c) => c.ok);
  return { ok, checks, browser };
}
