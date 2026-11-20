import { UsageError } from "./args.mjs";
import { request } from "./http.mjs";

const METHODS = ["GET", "POST", "PATCH", "PUT", "DELETE"];

export function apiPlan(cfg, method, urlPath, data) {
  const m = String(method ?? "").toUpperCase();
  if (!METHODS.includes(m)) throw new UsageError(`unknown method ${method}`, `use one of ${METHODS.join(", ")}, e.g. \`bin/verify api GET /api/samples\``);
  if (!urlPath || !urlPath.startsWith("/api/")) throw new UsageError(`path must start with /api/, got ${urlPath ?? "nothing"}`);
  let body;
  if (data !== undefined) {
    try {
      body = JSON.parse(data);
    } catch (e) {
      throw new UsageError(`--data is not valid JSON: ${e.message}`, `quote it for the shell: --data '{"status":"rejected"}'`);
    }
  }
  return { method: m, url: new URL(urlPath, cfg.apiBase).toString(), path: urlPath, body };
}

export async function api(cfg, method, urlPath, flags) {
  const plan = apiPlan(cfg, method, urlPath, flags.data);
  if (flags["dry-run"]) return { dryRun: true, ...plan };
  let res;
  try {
    res = await request(cfg.apiBase, plan.method, plan.path, { data: plan.body });
  } catch {
    throw new UsageError(`backend not answering at ${cfg.apiBase}`, "run `bin/verify doctor`, then `bin/verify seed`");
  }
  return { ok: res.ok, status: res.status, method: plan.method, path: plan.path, body: res.body };
}
