import { describe, it, expect, beforeEach, vi } from "vitest";

describe("getClientId", () => {
  beforeEach(() => sessionStorage.clear());

  it("mints a UUID on first call and persists to sessionStorage", async () => {
    const { getClientId } = await import("../src/lib/clientId");
    const id = getClientId();
    expect(id).toMatch(/^[0-9a-f-]{36}$/);
    expect(sessionStorage.getItem("himalaya.client_id")).toBe(id);
  });

  it("returns the same id on subsequent calls within a session", async () => {
    const { getClientId } = await import("../src/lib/clientId");
    expect(getClientId()).toBe(getClientId());
  });

  it("reuses an existing id in sessionStorage", async () => {
    sessionStorage.setItem("himalaya.client_id", "preset-uuid");
    vi.resetModules();
    const { getClientId } = await import("../src/lib/clientId");
    expect(getClientId()).toBe("preset-uuid");
  });

  it("falls back to a non-crypto UUID when crypto.randomUUID throws", async () => {
    // Simulates plain-http LAN deployment where randomUUID throws TypeError.
    const original = crypto.randomUUID;
    vi.stubGlobal("crypto", {
      ...crypto,
      randomUUID: () => {
        throw new TypeError("randomUUID requires a secure context");
      },
    });
    vi.resetModules();
    try {
      const { getClientId } = await import("../src/lib/clientId");
      const id = getClientId();
      expect(id).toMatch(/^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$/);
      expect(sessionStorage.getItem("himalaya.client_id")).toBe(id);
    } finally {
      vi.stubGlobal("crypto", { ...crypto, randomUUID: original });
    }
  });
});
