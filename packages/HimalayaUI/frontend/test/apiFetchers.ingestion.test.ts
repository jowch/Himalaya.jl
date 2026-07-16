// test/apiFetchers.ingestion.test.ts
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import * as api from "../src/api";

function mockFetchOnce(body: unknown, ok = true, status = 200): void {
  vi.spyOn(globalThis, "fetch").mockResolvedValueOnce({
    ok, status,
    text: async () => JSON.stringify(body),
    json: async () => body,
  } as Response);
}

describe("ingestion fetchers (Phase E1)", () => {
  beforeEach(() => vi.restoreAllMocks());
  afterEach(() => vi.restoreAllMocks());

  it("createExperiment POSTs {path,name?,patterns?} to /api/experiments", async () => {
    const spy = vi.spyOn(globalThis, "fetch").mockResolvedValueOnce({
      ok: true, status: 200,
      text: async () => JSON.stringify({ id: 7 }),
      json: async () => ({ id: 7 }),
    } as Response);
    const out = await api.createExperiment({ path: "/d", name: "X" }, { username: "u" });
    expect(out.id).toBe(7);
    const [url, init] = spy.mock.calls[0]!;
    expect(url).toBe("/api/experiments");
    expect(init!.method).toBe("POST");
    expect(JSON.parse(init!.body as string)).toEqual({ path: "/d", name: "X" });
    expect((init!.headers as Record<string, string>)["X-Username"]).toBe("u");
  });

  it("triggerScan POSTs to /api/experiments/:id/scan", async () => {
    const spy = mockFetchSpy({ id: 7, ingest_status: "scanning" });
    await api.triggerScan(7);
    expect(spy.mock.calls[0]![0]).toBe("/api/experiments/7/scan");
    expect((spy.mock.calls[0]![1] as RequestInit).method).toBe("POST");
  });

  it("listLoads GETs /api/experiments/:id/loads (nested roll-up)", async () => {
    const spy = mockFetchSpy([{
      load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null,
      frame_count: 0, note: null, samples: [],
    }]);
    const loads = await api.listLoads(7);
    expect(loads).toHaveLength(1);
    expect(loads[0]!.samples).toEqual([]);
    expect(spy.mock.calls[0]![0]).toBe("/api/experiments/7/loads");
  });

  it("updateExperiment PATCHes a geometry partial", async () => {
    const spy = mockFetchSpy({ id: 7, flight_path_m: 1.81, flight_path_m_source: "user" });
    const patch: api.ExperimentPatch = { flight_path_m: 1.81 };
    await api.updateExperiment(7, patch, { username: "u" });
    expect(spy.mock.calls[0]![0]).toBe("/api/experiments/7");
    expect((spy.mock.calls[0]![1] as RequestInit).method).toBe("PATCH");
    expect(JSON.parse((spy.mock.calls[0]![1] as RequestInit).body as string)).toEqual({ flight_path_m: 1.81 });
  });

  it("updateExperiment accepts name + description + file patterns", async () => {
    const spy = mockFetchSpy({ id: 7, name: "renamed" });
    const patch: api.ExperimentPatch = {
      name: "renamed", description: "AgBe run",
      image_pattern: "*.tif", metadata_pattern: "*.prp", integration_pattern: "*.dat",
    };
    await api.updateExperiment(7, patch, { username: "u" });
    const sent = JSON.parse((spy.mock.calls[0]![1] as RequestInit).body as string);
    expect(sent.name).toBe("renamed");
    expect(sent.image_pattern).toBe("*.tif");
  });

  it("suggestPaths GETs /api/fs/suggest?prefix=…", async () => {
    const spy = mockFetchSpy({ suggestions: ["/Volumes/data/ssrl/2026_04/1p7m"] });
    const out = await api.suggestPaths("/Volumes/data/ssrl/2026_04/1");
    expect(out.suggestions[0]).toContain("1p7m");
    expect(spy.mock.calls[0]![0]).toContain("/api/fs/suggest?");
    expect(spy.mock.calls[0]![0]).toContain("prefix=");
  });

  it("validatePath GETs /api/fs/validate?path=… and returns counts", async () => {
    const spy = mockFetchSpy({ ok: true, matched: 682, scanned: 700, message: null });
    const out = await api.validatePath("/Volumes/data/ssrl/2026_04/1p7m");
    expect(out.matched).toBe(682);
    expect(spy.mock.calls[0]![0]).toContain("/api/fs/validate?");
  });

  it("fetchManifest single-encodes path with slashes and spaces", async () => {
    const spy = mockFetchSpy({
      total: 0, matched: { image: 0, metadata: 0, integration: 0 }, unmatched: [],
    });
    await api.fetchManifest("/data/run 42");
    const url = spy.mock.calls[0]![0] as string;
    // URLSearchParams encodes '/' as %2F and ' ' as '+' (or %20).
    // Double-encoding would produce %252F for every slash — a clear signature.
    expect(url).toContain("/api/fs/manifest?");
    expect(url).not.toContain("%252F"); // no double-encoding
    // The path must be present in encoded form — either %2F or %2f for '/'
    expect(url.toLowerCase()).toContain("path=%2fdata%2frun");
  });

  it("fetchManifest single-encodes optional pattern params", async () => {
    const spy = mockFetchSpy({
      total: 0, matched: { image: 0, metadata: 0, integration: 0 }, unmatched: [],
    });
    await api.fetchManifest("/data/run 42", {
      image: "*.tif",
      metadata: "*.prp",
      integration: "*.dat",
    });
    const url = spy.mock.calls[0]![0] as string;
    expect(url).not.toContain("%252F");
    expect(url).toContain("image_pattern=");
    expect(url).toContain("metadata_pattern=");
    expect(url).toContain("integration_pattern=");
  });
});

function mockFetchSpy(body: unknown) {
  return vi.spyOn(globalThis, "fetch").mockResolvedValueOnce({
    ok: true, status: 200,
    text: async () => JSON.stringify(body),
    json: async () => body,
  } as Response);
}
