import fs from "node:fs";
import path from "node:path";
import { parseCsv } from "./csv.mjs";
import { getJson } from "./http.mjs";

export const CAPILLARY_DIR = "capillary";
export const AUTOSAMPLER_DIR = "autosampler";

export const capillaryDir = (cfg) => path.join(cfg.demoDir, CAPILLARY_DIR);
export const autosamplerDir = (cfg) => path.join(cfg.demoDir, AUTOSAMPLER_DIR);

export const STALE_LAYOUT_HINT = "run `bin/verify seed --reset` (or `make demo-clean demo-data`) to regenerate the demo in its current layout";

export function staleDemoLayout(cfg) {
  return fs.existsSync(path.join(cfg.demoDir, "TRUTH.csv")) &&
    fs.existsSync(path.join(cfg.demoDir, "data")) &&
    !fs.existsSync(capillaryDir(cfg));
}

export function readTruth(cfg, file = path.join(cfg.demoDir, "TRUTH.csv")) {
  if (!fs.existsSync(file)) return null;
  return parseCsv(fs.readFileSync(file, "utf8"));
}

export function readAutosamplerTruth(cfg, truthFile = path.join(cfg.demoDir, "TRUTH.csv")) {
  return readTruth(cfg, path.join(path.dirname(truthFile), "TRUTH_autosampler.csv"));
}

export async function findDemoExperiment(cfg) {
  const experiments = await getJson(cfg.apiBase, "/api/experiments");
  const dataDir = path.join(capillaryDir(cfg), "data");
  return experiments.find((e) => (e.setups ?? []).some((s) => s.data_dir && path.resolve(s.data_dir) === dataDir)) ?? null;
}
