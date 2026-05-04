import { describe, it, expect, beforeEach, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import {
  attachPersistence,
  rehydrate,
  STORAGE_KEY,
  SCHEMA_VERSION,
  type MutatorResolver,
} from "../../src/lib/queue/persistence";
import type { Mutator, OpKind } from "../../src/lib/queue/types";
import { makeFakeMutation } from "./helpers";

describe("persistence: mirrorToSessionStorage", () => {
  beforeEach(() => {
    sessionStorage.clear();
  });

  it("mirrors a pending mutation to sessionStorage on add", () => {
    const qc = new QueryClient();
    const mc = qc.getMutationCache();
    attachPersistence(mc);

    mc.add(makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_added", clientOpId: "op-1", payload: { q: 1.0 }, exposureId: 42 },
    }));

    const stored = JSON.parse(sessionStorage.getItem(STORAGE_KEY) ?? "[]");
    expect(stored).toHaveLength(1);
    expect(stored[0]).toMatchObject({
      schemaVersion: SCHEMA_VERSION,
      kind: "peak_added",
      clientOpId: "op-1",
    });
  });

  it("removes confirmed mutations from sessionStorage on update to success", () => {
    const qc = new QueryClient();
    const mc = qc.getMutationCache();
    attachPersistence(mc);

    const m = makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_added", clientOpId: "op-2", payload: { q: 2.0 } },
    });
    mc.add(m);
    expect(JSON.parse(sessionStorage.getItem(STORAGE_KEY) ?? "[]")).toHaveLength(1);

    // Transition to success and re-notify the cache.
    (m.state as { status: string }).status = "success";
    mc.notify({ type: "updated", mutation: m, action: { type: "success", data: {} } } as any);

    expect(JSON.parse(sessionStorage.getItem(STORAGE_KEY) ?? "[]")).toHaveLength(0);
  });

  it("attachPersistence returns an unsubscribe function that stops mirroring", () => {
    const qc = new QueryClient();
    const mc = qc.getMutationCache();
    const unsub = attachPersistence(mc);
    unsub();

    mc.add(makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_added", clientOpId: "op-3", payload: { q: 3.0 } },
    }));
    expect(JSON.parse(sessionStorage.getItem(STORAGE_KEY) ?? "[]")).toHaveLength(0);
  });
});

describe("persistence: rehydrate", () => {
  beforeEach(() => {
    sessionStorage.clear();
  });

  it("re-fires persisted ops via the matching mutator's request", async () => {
    const requestSpy = vi.fn().mockResolvedValue({ ok: true });
    const onMutate = vi.fn(() => ({ restore: () => {} }));
    const onSuccess = vi.fn();
    const mutator: Mutator<any, any, any> = {
      kind: "peak_added",
      onMutate,
      request: requestSpy,
      onSuccess,
    };
    sessionStorage.setItem(STORAGE_KEY, JSON.stringify([
      {
        schemaVersion: SCHEMA_VERSION,
        kind: "peak_added",
        clientOpId: "op-rehyd",
        payload: { q: 5.5 },
      },
    ]));

    const qc = new QueryClient();
    await rehydrate(qc, new Map<OpKind, Mutator<any, any, any>>([["peak_added", mutator]]));

    expect(onMutate).toHaveBeenCalledTimes(1);
    expect(requestSpy).toHaveBeenCalledTimes(1);
    expect(sessionStorage.getItem(STORAGE_KEY)).toBeNull();
  });

  it("drops persisted ops with mismatched schemaVersion and emits a toast count", async () => {
    sessionStorage.setItem(STORAGE_KEY, JSON.stringify([
      { schemaVersion: 999, kind: "peak_added", clientOpId: "op-bad", payload: {} },
      { schemaVersion: SCHEMA_VERSION, kind: "peak_added", clientOpId: "op-good", payload: { q: 1.0 } },
    ]));

    const requestSpy = vi.fn().mockResolvedValue({ ok: true });
    const mutator: Mutator<any, any, any> = {
      kind: "peak_added",
      onMutate: () => ({ restore: () => {} }),
      request: requestSpy,
      onSuccess: () => {},
    };
    const qc = new QueryClient();
    const result = await rehydrate(qc, new Map([["peak_added", mutator]]));

    expect(requestSpy).toHaveBeenCalledTimes(1); // only the good one
    expect(result.dropped).toBe(1);
    expect(result.replayed).toBe(1);
    expect(result.failed).toBe(0);
  });

  it("drops persisted ops whose kind has no registered mutator", async () => {
    sessionStorage.setItem(STORAGE_KEY, JSON.stringify([
      { schemaVersion: SCHEMA_VERSION, kind: "peak_added", clientOpId: "op-x", payload: {} },
    ]));
    const qc = new QueryClient();
    const result = await rehydrate(qc, new Map());  // empty mutator map
    expect(result.dropped).toBe(1);
    expect(result.replayed).toBe(0);
  });

  it("returns {dropped: 0, replayed: 0, failed: 0} when no persisted state", async () => {
    const qc = new QueryClient();
    const result = await rehydrate(qc, new Map());
    expect(result).toEqual({ dropped: 0, replayed: 0, failed: 0 });
  });

  it("clears the storage key on malformed JSON", async () => {
    sessionStorage.setItem(STORAGE_KEY, "{not json");
    const qc = new QueryClient();
    const result = await rehydrate(qc, new Map());
    expect(sessionStorage.getItem(STORAGE_KEY)).toBeNull();
    expect(result).toEqual({ dropped: 0, replayed: 0, failed: 0 });
  });

  it("rehydrate accepts a resolver function for dual-scope kinds", async () => {
    sessionStorage.setItem(
      STORAGE_KEY,
      JSON.stringify([
        {
          schemaVersion: SCHEMA_VERSION,
          kind: "add_tag",
          clientOpId: "op-resolver",
          payload: { experimentId: 1, sampleId: 10, key: "k", value: "v" },
        },
      ]),
    );
    const requestSpy = vi.fn().mockResolvedValue({ id: 5 });
    const sampleScopedMutator: Mutator<any, any, any> = {
      kind: "add_tag",
      onMutate: () => ({ restore: () => {} }),
      request: requestSpy,
      onSuccess: () => {},
    };
    const exposureScopedSpy = vi.fn();
    const exposureScopedMutator: Mutator<any, any, any> = {
      kind: "add_tag",
      onMutate: () => ({ restore: () => {} }),
      request: exposureScopedSpy,
      onSuccess: () => {},
    };
    const resolver: MutatorResolver = (op) => {
      if (op.kind !== "add_tag") return undefined;
      const p = op.payload as { experimentId?: number };
      return p?.experimentId !== undefined
        ? sampleScopedMutator
        : exposureScopedMutator;
    };
    const qc = new QueryClient();
    const result = await rehydrate(qc, resolver);
    expect(requestSpy).toHaveBeenCalledTimes(1);
    expect(exposureScopedSpy).not.toHaveBeenCalled();
    expect(result.replayed).toBe(1);
    expect(result.dropped).toBe(0);
  });

  it("rehydrate resolver receiving undefined returns dropped", async () => {
    sessionStorage.setItem(
      STORAGE_KEY,
      JSON.stringify([
        {
          schemaVersion: SCHEMA_VERSION,
          kind: "peak_added",
          clientOpId: "op-no-mut",
          payload: {},
        },
      ]),
    );
    const qc = new QueryClient();
    const result = await rehydrate(qc, () => undefined);
    expect(result.dropped).toBe(1);
    expect(result.replayed).toBe(0);
  });

  it("counts onMutate throws as failed without aborting remaining ops", async () => {
    sessionStorage.setItem(STORAGE_KEY, JSON.stringify([
      { schemaVersion: SCHEMA_VERSION, kind: "peak_added", clientOpId: "op-throw", payload: { bad: true } },
      { schemaVersion: SCHEMA_VERSION, kind: "peak_added", clientOpId: "op-good", payload: { q: 1.0 } },
    ]));
    const requestSpy = vi.fn().mockResolvedValue({ ok: true });
    const onMutate = vi.fn((p: { bad?: boolean }) => {
      if (p?.bad) throw new Error("malformed payload");
      return { restore: () => {} };
    });
    const mutator: Mutator<any, any, any> = {
      kind: "peak_added",
      onMutate,
      request: requestSpy,
      onSuccess: () => {},
    };
    const qc = new QueryClient();
    const result = await rehydrate(qc, new Map([["peak_added", mutator]]));
    // First op's onMutate threw → counted as failed; loop continued so second op still ran.
    expect(result.failed).toBe(1);
    expect(result.replayed).toBe(1);
    expect(requestSpy).toHaveBeenCalledTimes(1);
  });

  it("invokes restore() when a replayed request rejects", async () => {
    const restore = vi.fn();
    sessionStorage.setItem(STORAGE_KEY, JSON.stringify([
      { schemaVersion: SCHEMA_VERSION, kind: "peak_added", clientOpId: "op-rollback", payload: {} },
    ]));
    const mutator: Mutator<any, any, any> = {
      kind: "peak_added",
      onMutate: () => ({ restore }),
      request: vi.fn().mockRejectedValue(new Error("boom")),
      onSuccess: () => {},
    };
    const qc = new QueryClient();
    const result = await rehydrate(qc, new Map([["peak_added", mutator]]));
    expect(restore).toHaveBeenCalledTimes(1);
    expect(result.failed).toBe(1);
  });

  it("counts failed replays separately from successful ones", async () => {
    sessionStorage.setItem(STORAGE_KEY, JSON.stringify([
      { schemaVersion: SCHEMA_VERSION, kind: "peak_added", clientOpId: "op-fail", payload: {} },
    ]));
    const mutator: Mutator<any, any, any> = {
      kind: "peak_added",
      onMutate: () => ({ restore: () => {} }),
      request: vi.fn().mockRejectedValue(new Error("boom")),
      onSuccess: () => {},
    };
    const qc = new QueryClient();
    const result = await rehydrate(qc, new Map([["peak_added", mutator]]));
    expect(result.failed).toBe(1);
    expect(result.replayed).toBe(0);
    expect(result.dropped).toBe(0);
  });
});
