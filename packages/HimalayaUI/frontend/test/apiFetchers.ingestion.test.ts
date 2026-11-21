// test/apiFetchers.ingestion.test.ts
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import * as api from "../src/api";

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

  it("suggestPaths GETs /api/fs/suggest?prefix=…", async () => {
    const spy = mockFetchSpy({ suggestions: ["/Volumes/data/ssrl/2026_04/1p7m"] });
    const out = await api.suggestPaths("/Volumes/data/ssrl/2026_04/1");
    expect(out.suggestions[0]).toContain("1p7m");
    expect(spy.mock.calls[0]![0]).toContain("/api/fs/suggest?");
    expect(spy.mock.calls[0]![0]).toContain("prefix=");
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
