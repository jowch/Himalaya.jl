import fs from "node:fs";
import path from "node:path";
import { UsageError } from "./args.mjs";

const SUFFIX = "_tot.dat";

export function lateDir(cfg) {
  return path.join(cfg.demoDir, "late");
}

function heldFiles(cfg) {
  const dir = lateDir(cfg);
  if (!fs.existsSync(dir)) return [];
  return fs.readdirSync(dir).filter((f) => f.endsWith(SUFFIX)).sort();
}

export function heldTraces(cfg) {
  return heldFiles(cfg).map((f) => f.slice(0, -SUFFIX.length));
}

export function deliverLate(cfg) {
  const files = heldFiles(cfg);
  if (!files.length) {
    throw new UsageError(`no held traces in ${lateDir(cfg)}`, "generate a demo with one: `bin/verify seed --reset --late-trace`");
  }
  const tot = path.join(cfg.demoDir, "analysis", "automatic_analysis", "tot_files");
  const now = new Date();
  for (const f of files) {
    const dst = path.join(tot, f);
    fs.renameSync(path.join(lateDir(cfg), f), dst);
    fs.utimesSync(dst, now, now);
  }
  fs.rmSync(lateDir(cfg), { recursive: true, force: true });
  return {
    delivered: files,
    next: "open the experiment (or run Rescan from the ⌘K palette) so a pass picks them up, then run check-truth",
  };
}
