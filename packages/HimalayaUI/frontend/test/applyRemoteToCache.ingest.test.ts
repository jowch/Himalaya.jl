// test/applyRemoteToCache.ingest.test.ts
import { describe, it, expect, vi } from "vitest";
import { QueryClient, QueryObserver } from "@tanstack/react-query";
import { applyRemoteToCache, invalidateIngestFrameCache } from "../src/lib/queue/applyRemoteToCache";
import { queryKeys } from "../src/queries";
import type { SseEvent } from "../src/lib/queue/types";

function ingestFrame(kind: string, experiment_id: number, extra: Record<string, unknown> = {}): SseEvent {
  // broadcast_progress! (events.jl:1158-1186) emits kind + payload at top-level;
  // it does NOT emit a top-level entity_id. The payload sub-object carries
  // experiment_id, processed, total. entity_id is -1 sentinel here (not emitted
  // by the real backend for ingest frames; the arms read payload.experiment_id).
  return {
    kind,
    entity_id: -1,
    entity_type: "experiment",
    payload: { experiment_id, ...extra },
  } as unknown as SseEvent;
}

describe("ingest_* cache arms (Phase E1)", () => {
  it("ingest_progress invalidates the experiment samples and the workspace sample list, never peaks/indices", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(ingestFrame("ingest_progress", 7, { processed: 100, total: 680 }), qc);
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.experimentCorpusSamples(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.samples(7)));
    // the default arm's peaks/indices invalidation must NOT fire
    expect(keys).not.toContain(JSON.stringify(queryKeys.peaks(7)));
    expect(keys).not.toContain(JSON.stringify(queryKeys.indices(7)));
  });

  it("ingest_complete invalidates the sample lists, the experiment detail and its problems", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(ingestFrame("ingest_complete", 7), qc);
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.experiment(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.problems(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.experimentCorpusSamples(7)));
  });

  it("ingest_complete refetches the experiment's labels, so role suggestions made at ingest show up", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(ingestFrame("ingest_complete", 7), qc);
    invalidateIngestFrameCache(qc, 8, false, "persist");
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.labels(7)));
    expect(keys).not.toContain(JSON.stringify(queryKeys.labels(8)));
  });

  it("an analyze_run whose trace changed refetches that exposure's trace and peaks, so a late trace shows its markers", async () => {
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    const calls = { trace: 0, peaks: 0, otherPeaks: 0 };
    const watch = (key: readonly unknown[], fn: () => Promise<unknown>) => new QueryObserver(qc, { queryKey: key, queryFn: fn }).subscribe(() => {});
    const stop = [
      watch(queryKeys.trace(11), () => { calls.trace += 1; return calls.trace === 1 ? Promise.reject(new Error("404")) : Promise.resolve({ q: [1], I: [1], sigma: [1] }); }),
      watch(queryKeys.peaks(11), () => { calls.peaks += 1; return Promise.resolve([]); }),
      watch(queryKeys.peaks(12), () => { calls.otherPeaks += 1; return Promise.resolve([]); }),
    ];
    await vi.waitFor(() => expect(calls).toEqual({ trace: 1, peaks: 1, otherPeaks: 1 }));
    await vi.waitFor(() => expect([queryKeys.trace(11), queryKeys.peaks(11), queryKeys.peaks(12)].map((k) => qc.getQueryState(k)?.fetchStatus)).toEqual(["idle", "idle", "idle"]));
    const run = (before: string | null, after: string) => ({
      kind: "analyze_run", entity_id: 11, entity_type: "exposure",
      payload: { trace_hash_before: before, trace_hash_after: after },
    }) as unknown as SseEvent;
    applyRemoteToCache(run("h1", "h1"), qc);
    await new Promise((r) => setTimeout(r, 20));
    expect(calls).toEqual({ trace: 1, peaks: 1, otherPeaks: 1 });
    applyRemoteToCache(run(null, "h1"), qc);
    await vi.waitFor(() => expect(calls).toEqual({ trace: 2, peaks: 2, otherPeaks: 1 }));
    await vi.waitFor(() => expect(qc.getQueryState(queryKeys.trace(11))?.status).toBe("success"));
    stop.forEach((s) => s());
  });

  it("ingest_failed invalidates the experiment detail and its problems, so Details shows the error", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(ingestFrame("ingest_failed", 7), qc);
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.experiment(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.problems(7)));
  });

  it("ingest_started / ingest_failed do not throw and do not poison peaks/indices", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(ingestFrame("ingest_started", 7), qc);
    applyRemoteToCache(ingestFrame("ingest_failed", 7), qc);
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).not.toContain(JSON.stringify(queryKeys.peaks(7)));
  });
});

describe("ingest_progress stage gating", () => {
  it("SKIPS the sample refetch during discovery", () => {
    // Nothing is committed during discovery — scan_and_group! inserts every
    // load/sample/exposure in ONE transaction that has not opened yet — so these
    // refetches can only return the same empty result, once per progress frame,
    // against the slowest payload on the page.
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(
      ingestFrame("ingest_progress", 7, { processed: 300, total: 1100, stage: "discovery" }),
      qc,
    );
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).not.toContain(JSON.stringify(queryKeys.experimentCorpusSamples(7)));
    expect(keys).not.toContain(JSON.stringify(queryKeys.samples(7)));
  });

  it("STILL refetches during analyzing (the persist txn has committed)", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(
      ingestFrame("ingest_progress", 7, { processed: 92, total: 604, stage: "analyzing" }),
      qc,
    );
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.experimentCorpusSamples(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.samples(7)));
  });

  it("keeps the original behavior when no stage is present", () => {
    // An older backend, or an ingest_started frame: must not silently lose the
    // refresh just because the field is missing.
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(ingestFrame("ingest_progress", 7, { processed: 1, total: 2 }), qc);
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.experimentCorpusSamples(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.samples(7)));
  });

  it("ingest_complete always refetches, even carrying the skipped stage", () => {
    // Called on the HELPER, not through applyRemoteToCache: that dispatcher's
    // ingest_complete arm passes no stage at all, so routing through it could not
    // exercise the isComplete-vs-stage precedence this asserts. App.tsx DOES
    // forward stage for all four ingest kinds, so the conflicting combination is
    // reachable in production.
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    invalidateIngestFrameCache(qc, 7, true, "discovery");
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.experimentCorpusSamples(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.samples(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.experiment(7)));
  });

  it("the discovery skip applies only to NON-terminal frames", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    invalidateIngestFrameCache(qc, 7, false, "discovery");
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).not.toContain(JSON.stringify(queryKeys.experimentCorpusSamples(7)));
    expect(keys).not.toContain(JSON.stringify(queryKeys.samples(7)));
  });
});

describe("selection cache arms", () => {
  const frame = (kind: string, entity_type: string, entity_id: number, payload: Record<string, unknown> = {}) =>
    ({ id: 1, kind, entity_type, entity_id, payload }) as unknown as SseEvent;
  const keysOf = (spy: { mock: { calls: unknown[][] } }) => spy.mock.calls.map((c) => JSON.stringify((c[0] as { queryKey?: unknown } | undefined)?.queryKey));

  it("an analyze_run on an average refetches its sample's images, so a new average replaces the waiting notice", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(frame("analyze_run", "exposure", 41, { trace_hash_before: null, trace_hash_after: "h", average_of_sample: 6 }), qc);
    expect(keysOf(spy)).toContain(JSON.stringify(["sample", 6, "exposures"]));
  });

  it("an analyze_run on a single image leaves the image lists alone", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(frame("analyze_run", "exposure", 13, { trace_hash_before: "a", trace_hash_after: "a", findpeaks_skipped: true }), qc);
    expect(keysOf(spy)).not.toContain(JSON.stringify(["sample", 6, "exposures"]));
  });

  it("a set_exposure_status leaves other samples' image lists invalidated, so a sample whose average changed still refetches", () => {
    const qc = new QueryClient();
    qc.setQueryData(["sample", 1, "exposures"], [{ id: 3, status: null }]);
    qc.setQueryData(["sample", 2, "exposures"], [{ id: 5, status: null }]);
    void qc.invalidateQueries({ queryKey: ["sample", 1, "exposures"] });
    applyRemoteToCache(frame("set_exposure_status", "exposure", 5, { status: "rejected" }), qc);
    expect(qc.getQueryState(["sample", 1, "exposures"])?.isInvalidated).toBe(true);
    expect(qc.getQueryData<{ status: string | null }[]>(["sample", 2, "exposures"])?.[0]?.status).toBe("rejected");
  });

  it("an update_sample leaves other experiments' sample lists untouched", () => {
    const qc = new QueryClient();
    qc.setQueryData(["experiment", 1, "samples"], [{ id: 1, name: "a" }]);
    qc.setQueryData(["experiment", 2, "samples"], [{ id: 9, name: "b" }]);
    void qc.invalidateQueries({ queryKey: ["experiment", 1, "samples"] });
    applyRemoteToCache(frame("update_sample", "sample", 9, { name: "c" }), qc);
    expect(qc.getQueryState(["experiment", 1, "samples"])?.isInvalidated).toBe(true);
    expect(qc.getQueryData<{ name: string }[]>(["experiment", 2, "samples"])?.[0]?.name).toBe("c");
  });

  it("a set_selection from another client refetches that sample's images and the sample lists", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(frame("set_selection", "sample", 6), qc);
    const keys = keysOf(spy);
    expect(keys).toContain(JSON.stringify(["sample", 6, "exposures"]));
    expect(keys).toContain(JSON.stringify(queryKeys.corpusSamples));
  });
});
