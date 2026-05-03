// Per-mutation idempotency token. Minted fresh at each mutate-call site
// (in contrast to clientId, which is per-tab and persisted in sessionStorage).
// Sent on every mutation as X-Client-Op-Id; backend uses it for request-level
// retry idempotency (Stripe-style) — same UUID returns the cached response,
// new UUID re-executes.

function mintUuid(): string {
  if (typeof crypto !== "undefined" && typeof crypto.randomUUID === "function") {
    try {
      return crypto.randomUUID();
    } catch {
      // fall through (plain-http deployments throw TypeError)
    }
  }
  return "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".replace(/[xy]/g, (c) => {
    const r = (Math.random() * 16) | 0;
    const v = c === "x" ? r : (r & 0x3) | 0x8;
    return v.toString(16);
  });
}

export function newClientOpId(): string {
  return mintUuid();
}
