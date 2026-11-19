import { spawnSync } from "node:child_process";

export const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

export function commandOf(pid) {
  const ps = spawnSync("ps", ["-o", "command=", "-p", String(pid)], { encoding: "utf8" });
  return ps.status === 0 ? ps.stdout.trim() : null;
}

export function alive(pid) {
  try {
    process.kill(pid, 0);
    return true;
  } catch {
    return false;
  }
}

export async function waitGone(pid, ms) {
  for (let t = 0; t < ms && alive(pid); t += 100) await sleep(100);
  return !alive(pid);
}
