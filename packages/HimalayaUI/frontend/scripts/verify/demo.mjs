import fs from "node:fs";
import path from "node:path";
import { parseCsv } from "./csv.mjs";
import { getJson } from "./http.mjs";

export function readTruth(cfg, file = path.join(cfg.demoDir, "TRUTH.csv")) {
  if (!fs.existsSync(file)) return null;
  return parseCsv(fs.readFileSync(file, "utf8"));
}

export async function findDemoExperiment(cfg) {
  const experiments = await getJson(cfg.apiBase, "/api/experiments");
  const dataDir = path.join(cfg.demoDir, "data");
  return experiments.find((e) => e.data_dir && path.resolve(e.data_dir) === dataDir) ?? null;
}
