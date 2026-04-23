import { describe, it, expect, vi, beforeEach } from "vitest";
import * as api from "../src/api";
import { ApiError } from "../src/api";

describe("api", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("sends X-Username header on mutating calls when opts.username provided", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ id: 1, username: "alice" }), { status: 200 }),
    );
    await api.createUser("alice", { username: "alice" });
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
  });

  it("omits X-Username on GET requests", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response("[]", { status: 200 }),
    );
    await api.listUsers();
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    const headers = (init.headers as Record<string, string>) ?? {};
    expect(headers["X-Username"]).toBeUndefined();
  });

  it("does not export setUsername/getUsername", () => {
    expect((api as unknown as Record<string, unknown>).setUsername).toBeUndefined();
    expect((api as unknown as Record<string, unknown>).getUsername).toBeUndefined();
  });

  it("throws ApiError on non-2xx responses with server error message", async () => {
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ error: "nope" }), { status: 400 }),
    );
    await expect(api.createUser("bob")).rejects.toMatchObject({
      name: "ApiError",
      status: 400,
      message: "nope",
    });
  });

  it("getTrace returns parsed q/I/sigma arrays", async () => {
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ q: [0.1, 0.2], I: [10, 20], sigma: [1, 2] }),
        { status: 200 }),
    );
    const t = await api.getTrace(42);
    expect(t.q).toEqual([0.1, 0.2]);
    expect(t.I).toEqual([10, 20]);
    expect(t.sigma).toEqual([1, 2]);
  });

  it("addPeak posts {q} with X-Username and returns parsed peak", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({
        id: 7, exposure_id: 42, q: 0.15, source: "manual", stale_indices: 3,
      }), { status: 201 }),
    );
    const p = await api.addPeak(42, 0.15, { username: "alice" });
    expect(p.id).toBe(7);
    expect(p.source).toBe("manual");
    expect(p.stale_indices).toBe(3);
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
    expect(init.body).toBe(JSON.stringify({ q: 0.15 }));
  });

  it("removePeak sends DELETE with X-Username", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(null, { status: 204 }),
    );
    await api.removePeak(7, { username: "alice" });
    const [url, init] = fetchSpy.mock.calls[0]! as [string, RequestInit];
    expect(url).toBe("/api/peaks/7");
    expect(init.method).toBe("DELETE");
    expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
  });

  it("reanalyzeExposure posts empty body and returns {id, analyzed}", async () => {
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ id: 42, analyzed: true }), { status: 200 }),
    );
    const r = await api.reanalyzeExposure(42, { username: "alice" });
    expect(r.id).toBe(42);
    expect(r.analyzed).toBe(true);
  });
});
