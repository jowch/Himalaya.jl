export async function request(base, method, urlPath, { data, timeoutMs = 10000 } = {}) {
  const init = { method, signal: AbortSignal.timeout(timeoutMs), headers: {} };
  if (data !== undefined) {
    init.headers["Content-Type"] = "application/json";
    init.body = typeof data === "string" ? data : JSON.stringify(data);
  }
  const res = await fetch(new URL(urlPath, base), init);
  const text = await res.text();
  let body = text;
  try {
    body = text ? JSON.parse(text) : null;
  } catch {
    body = text;
  }
  return { status: res.status, ok: res.ok, body };
}

export async function getJson(base, urlPath, opts) {
  const r = await request(base, "GET", urlPath, opts);
  if (!r.ok) throw new Error(`GET ${urlPath} returned ${r.status}`);
  return r.body;
}

export async function reachable(base, urlPath, timeoutMs = 2000) {
  try {
    return await request(base, "GET", urlPath, { timeoutMs });
  } catch {
    return null;
  }
}
