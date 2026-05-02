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
});
