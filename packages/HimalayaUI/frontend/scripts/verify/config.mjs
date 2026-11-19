import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
export const REPO = path.resolve(here, "../../../../..");
export const FRONTEND = path.join(REPO, "packages/HimalayaUI/frontend");

export function loadConfig(env = process.env) {
  const demoDir = path.resolve(REPO, env.DEMO_DIR ?? "scratch/demo-experiment");
  const demoDbDir = path.resolve(REPO, env.DEMO_DB_DIR ?? "scratch/demo");
  const stateDir = path.resolve(REPO, env.VERIFY_STATE_DIR ?? "scratch/verify");
  const apiPort = Number(env.DEMO_API_PORT ?? 8080);
  const uiPort = Number(env.DEMO_UI_PORT ?? 5173);
  const cdpPort = Number(env.VERIFY_CDP_PORT ?? 9333);
  return {
    demoDir,
    demoDbDir,
    stateDir,
    apiPort,
    uiPort,
    cdpPort,
    apiBase: `http://127.0.0.1:${apiPort}`,
    uiBase: `http://127.0.0.1:${uiPort}`,
    cdpBase: `http://127.0.0.1:${cdpPort}`,
  };
}
