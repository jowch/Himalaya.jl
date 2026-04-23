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

  it("listIndices returns indices with predicted_q and enriched peaks", async () => {
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify([
        {
          id: 1, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 1.0,
          r_squared: 0.99, lattice_d: 12.5, status: "candidate",
          predicted_q: [0.7071, 0.866, 1.0],
          peaks: [{ peak_id: 10, ratio_position: 1, residual: 0.001, q_observed: 0.71 }],
        },
      ]), { status: 200 }),
    );
    const indices = await api.listIndices(42);
    expect(indices).toHaveLength(1);
    expect(indices[0]!.predicted_q).toEqual([0.7071, 0.866, 1.0]);
    expect(indices[0]!.peaks[0]!.q_observed).toBeCloseTo(0.71);
  });

  it("listGroups fetches groups for exposure", async () => {
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify([
        { id: 1, exposure_id: 42, kind: "auto",   active: false, members: [10] },
        { id: 2, exposure_id: 42, kind: "custom", active: true,  members: [10, 11] },
      ]), { status: 200 }),
    );
    const groups = await api.listGroups(42);
    expect(groups).toHaveLength(2);
    expect(groups[1]!.kind).toBe("custom");
    expect(groups[1]!.active).toBe(true);
  });

  it("addIndexToGroup posts {index_id} with X-Username", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({
        id: 2, exposure_id: 42, kind: "custom", active: true, members: [10, 11],
      }), { status: 200 }),
    );
    const g = await api.addIndexToGroup(2, 11, { username: "alice" });
    expect(g.members).toEqual([10, 11]);
    const [, init] = fetchSpy.mock.calls[0]! as [string, RequestInit];
    expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
    expect(init.body).toBe(JSON.stringify({ index_id: 11 }));
  });

  it("removeIndexFromGroup sends DELETE with X-Username", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({
        id: 2, exposure_id: 42, kind: "custom", active: true, members: [10],
      }), { status: 200 }),
    );
    await api.removeIndexFromGroup(2, 11, { username: "alice" });
    const [url, init] = fetchSpy.mock.calls[0]! as [string, RequestInit];
    expect(url).toBe("/api/groups/2/members/11");
    expect(init.method).toBe("DELETE");
    expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
  });
});
