import { UsageError } from "./args.mjs";

export function parseThrottle(value) {
  if (value === undefined || value.trim() === "") return null;
  const m = /^\s*(\d+(?:\.\d+)?)\s*mbit\s*,\s*(\d+)\s*ms\s*$/i.exec(value);
  if (!m || Number(m[1]) <= 0) {
    throw new UsageError(`VERIFY_THROTTLE must look like 20mbit,30ms, got ${value}`, "set it to the link you want to emulate, for example VERIFY_THROTTLE=20mbit,30ms for an SSH tunnel from a laptop");
  }
  return { mbit: Number(m[1]), latencyMs: Number(m[2]) };
}

export async function applyThrottle(context, page, throttle) {
  if (!throttle) return;
  const cdp = await context.newCDPSession(page);
  const bytesPerSecond = (throttle.mbit * 1e6) / 8;
  await cdp.send("Network.enable");
  await cdp.send("Network.emulateNetworkConditions", {
    offline: false,
    latency: throttle.latencyMs,
    downloadThroughput: bytesPerSecond,
    uploadThroughput: bytesPerSecond,
  });
}
