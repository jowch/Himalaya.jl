import { describe, it, expect } from "vitest";
import { newClientOpId } from "../src/lib/clientOpId";

describe("newClientOpId", () => {
  it("mints a UUID v4-shaped string", () => {
    const id = newClientOpId();
    expect(id).toMatch(/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/);
  });

  it("returns unique IDs across calls", () => {
    const ids = new Set<string>();
    for (let i = 0; i < 100; i++) ids.add(newClientOpId());
    expect(ids.size).toBe(100);
  });

  it("does NOT persist to storage (per-call, not per-tab)", () => {
    const before = sessionStorage.length;
    newClientOpId();
    expect(sessionStorage.length).toBe(before);
  });
});
