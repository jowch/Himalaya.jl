// test/apiFetchers.ingestion.test.ts
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import * as api from "../src/api";

describe("ingestion fetchers (Phase E1)", () => {
  beforeEach(() => vi.restoreAllMocks());
  afterEach(() => vi.restoreAllMocks());

  it("createExperiment POSTs {path,name?,setups} to /api/experiments", async () => {
    const spy = vi.spyOn(globalThis, "fetch").mockResolvedValueOnce({
      ok: true, status: 200,
      text: async () => JSON.stringify({ id: 7 }),
      json: async () => ({ id: 7 }),
    } as Response);
    const body = { path: "/d", name: "X", setups: [{ root: "/d/setup1", data_dir: "/d/setup1/data", included: true }] };
    const out = await api.createExperiment(body, { username: "u" });
    expect(out.id).toBe(7);
    const [url, init] = spy.mock.calls[0]!;
    expect(url).toBe("/api/experiments");
    expect(init!.method).toBe("POST");
    expect(JSON.parse(init!.body as string)).toEqual(body);
    expect((init!.headers as Record<string, string>)["X-Username"]).toBe("u");
  });

  it("triggerScan POSTs to /api/experiments/:id/scan", async () => {
    const spy = mockFetchSpy({ id: 7, ingest_status: "scanning" });
    await api.triggerScan(7);
    expect(spy.mock.calls[0]![0]).toBe("/api/experiments/7/scan");
    expect((spy.mock.calls[0]![1] as RequestInit).method).toBe("POST");
  });

  it("listFolder GETs /api/fs/list?path=… with the path encoded once", async () => {
    const spy = mockFetchSpy({ path: "/data/run 42", parent: "/data", folders: ["setup1"] });
    const out = await api.listFolder("/data/run 42");
    expect(out.folders).toEqual(["setup1"]);
    expect(spy.mock.calls[0]![0]).toBe("/api/fs/list?path=%2Fdata%2Frun%2042");
  });

  it("describeFolder GETs /api/fs/describe?path=…", async () => {
    const spy = mockFetchSpy({ path: "/data/beamtime-01", kind: "beamtime", setups: 2, added: null });
    const out = await api.describeFolder("/data/beamtime-01");
    expect(out.kind).toBe("beamtime");
    expect(spy.mock.calls[0]![0]).toBe("/api/fs/describe?path=%2Fdata%2Fbeamtime-01");
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
