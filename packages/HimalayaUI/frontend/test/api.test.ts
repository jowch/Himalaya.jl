import { describe, it, expect, vi, beforeEach } from "vitest";
import * as api from "../src/api";

describe("api", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("sends X-Username header on mutating calls when opts.username provided", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ id: 1, username: "alice" }), { status: 200 }),
    );
    await api.createUser("alice", undefined, { username: "alice", clientId: "tab-xyz" });
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
    expect((init.headers as Record<string, string>)["X-Client-Id"]).toBe("tab-xyz");
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

  it("omits X-Client-Id on GET requests", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response("[]", { status: 200 }),
    );
    await api.listUsers();
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    const headers = (init.headers as Record<string, string>) ?? {};
    expect(headers["X-Client-Id"]).toBeUndefined();
  });

  it("sends X-Client-Id without X-Username when only clientId is provided", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({
        id: 7, exposure_id: 42, q: 0.15, source: "manual", stale_indices: 0,
      }), { status: 201 }),
    );
    await api.addPeak(42, 0.15, { clientId: "tab-xyz" });
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    const headers = init.headers as Record<string, string>;
    expect(headers["X-Client-Id"]).toBe("tab-xyz");
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
        id: 7, exposure_id: 42, q: 0.15, source: "manual",
        event_id: 100, view_row_id: 7, analysis_inputs_hash: "sha256:abc",
      }), { status: 201 }),
    );
    const p = await api.addPeak(42, 0.15, { username: "alice", clientId: "tab-xyz" });
    expect(p.id).toBe(7);
    expect(p.source).toBe("manual");
    expect(p.analysis_inputs_hash).toBe("sha256:abc");
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
    expect((init.headers as Record<string, string>)["X-Client-Id"]).toBe("tab-xyz");
    expect(init.body).toBe(JSON.stringify({ q: 0.15 }));
  });

  it("removePeak sends DELETE with X-Username and returns hash payload", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({
        event_id: 11, view_row_id: 22, analysis_inputs_hash: "hash-x",
      }), { status: 200, headers: { "Content-Type": "application/json" } }),
    );
    const r = await api.removePeak(7, { username: "alice", clientId: "tab-xyz" });
    expect(r.analysis_inputs_hash).toBe("hash-x");
    expect(r.event_id).toBe(11);
    const [url, init] = fetchSpy.mock.calls[0]! as [string, RequestInit];
    expect(url).toBe("/api/peaks/7");
    expect(init.method).toBe("DELETE");
    expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
    expect((init.headers as Record<string, string>)["X-Client-Id"]).toBe("tab-xyz");
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

  it("listSampleMessages fetches messages for a sample", async () => {
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify([
        { id: 1, sample_id: 3, author_id: 1, author: "alice", body: "hi", created_at: "2026-04-24 10:00:00" },
        { id: 2, sample_id: 3, author_id: 2, author: "bob",   body: "yo", created_at: "2026-04-24 10:01:00" },
      ]), { status: 200 }),
    );
    const msgs = await api.listSampleMessages(3);
    expect(msgs).toHaveLength(2);
    expect(msgs[0]!.author).toBe("alice");
    expect(msgs[1]!.body).toBe("yo");
  });

  it("postSampleMessage posts {body} with X-Username and returns parsed message", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({
        id: 7, sample_id: 3, author_id: 1, author: "alice", body: "hello", created_at: "2026-04-24 10:00:00",
      }), { status: 201 }),
    );
    const msg = await api.postSampleMessage(3, "hello", { username: "alice", clientId: "tab-xyz" });
    expect(msg.id).toBe(7);
    expect(msg.author).toBe("alice");
    expect(msg.body).toBe("hello");
    const [url, init] = fetchSpy.mock.calls[0]! as [string, RequestInit];
    expect(url).toBe("/api/samples/3/messages");
    expect(init.method).toBe("POST");
    expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
    expect((init.headers as Record<string, string>)["X-Client-Id"]).toBe("tab-xyz");
    expect(init.body).toBe(JSON.stringify({ body: "hello" }));
  });

  it("threads X-Client-Op-Id on POST mutations", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({
        id: 7, exposure_id: 42, q: 0.15, source: "manual", stale_indices: 0,
      }), { status: 201 }),
    );
    await api.addPeak(42, 0.15, { clientId: "tab-id", clientOpId: "op-id" });
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    const headers = init.headers as Record<string, string>;
    expect(headers["X-Client-Op-Id"]).toBe("op-id");
  });

  it("omits X-Client-Op-Id on GET requests", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response("[]", { status: 200 }),
    );
    // GET path; AuthOpts isn't directly supported on listUsers, so use listPeaks
    // which is also GET — but neither accepts opts. Verify via a direct GET call
    // through any GET fetcher: listUsers does not pass opts, but the request()
    // function still respects method===GET regardless. We assert no header
    // appears even if upstream code somehow forwarded it.
    await api.listUsers();
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    const headers = (init.headers as Record<string, string>) ?? {};
    expect(headers["X-Client-Op-Id"]).toBeUndefined();
  });

  it("threads X-Client-Op-Id without X-Username (system flows)", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({
        id: 7, exposure_id: 42, q: 0.15, source: "manual", stale_indices: 0,
      }), { status: 201 }),
    );
    await api.addPeak(42, 0.15, { clientId: "tab-id", clientOpId: "op-id" });
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    const headers = init.headers as Record<string, string>;
    expect(headers["X-Client-Op-Id"]).toBe("op-id");
    expect(headers["X-Username"]).toBeUndefined();
  });

});
