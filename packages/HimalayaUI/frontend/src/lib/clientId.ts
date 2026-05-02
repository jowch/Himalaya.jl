const KEY = "himalaya.client_id";

// `crypto.randomUUID` is only defined in secure contexts (https:// or
// http://localhost). HimalayaUI also targets plain-http LAN deployments
// (e.g. http://lab-host.local:8080), where calling it throws TypeError. Since
// `getClientId` runs at module-load time via `queries.ts`, that throw would
// crash the SPA before mount. Fall back to a non-crypto UUID v4 — collision
// risk is irrelevant because the id is only ever compared within a single
// browser's set of open tabs.
function mintUuid(): string {
  if (typeof crypto !== "undefined" && typeof crypto.randomUUID === "function") {
    try {
      return crypto.randomUUID();
    } catch {
      // fall through to non-crypto fallback
    }
  }
  return "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".replace(/[xy]/g, (c) => {
    const r = (Math.random() * 16) | 0;
    const v = c === "x" ? r : (r & 0x3) | 0x8;
    return v.toString(16);
  });
}

export function getClientId(): string {
  let id = sessionStorage.getItem(KEY);
  if (!id) {
    id = mintUuid();
    sessionStorage.setItem(KEY, id);
  }
  return id;
}
