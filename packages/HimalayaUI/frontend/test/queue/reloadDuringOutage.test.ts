import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import {
  attachPersistence,
  rehydrate,
  STORAGE_KEY,
  SCHEMA_VERSION,
} from "../../src/lib/queue/persistence";
import type { Mutator, OpKind } from "../../src/lib/queue/types";
import { pendingDeferreds } from "../../src/lib/queue/deferred";
import { NetworkError } from "../../src/api";
import * as toast from "../../src/lib/toast";

function storeOps(...ids: string[]) {
  sessionStorage.setItem(
    STORAGE_KEY,
    JSON.stringify(
      ids.map((id) => ({
        schemaVersion: SCHEMA_VERSION,
        kind: "set_exposure_status",
        clientOpId: id,
        payload: { kind: "set_exposure_status", clientOpId: id, exposureId: 7, status: "rejected" },
      })),
    ),
  );
}

function stored(): { clientOpId: string }[] {
  return JSON.parse(sessionStorage.getItem(STORAGE_KEY) ?? "[]");
}

function downFor(ms: number) {
  const upAt = Date.now() + ms;
  return vi.fn(async (_p: { clientOpId: string }) => {
    if (Date.now() < upAt) throw new NetworkError("PUT", "/api/x");
    return { ok: true };
  });
}

describe("an edit queued when the page reloads during an outage", () => {
  let qc: QueryClient;
  let toastSpy: ReturnType<typeof vi.spyOn>;

  beforeEach(() => {
    vi.useFakeTimers();
    sessionStorage.clear();
    pendingDeferreds.clear();
    qc = new QueryClient();
    toastSpy = vi.spyOn(toast, "showToast").mockImplementation(() => {});
  });
  afterEach(() => {
    toastSpy.mockRestore();
    vi.useRealTimers();
  });

  function mutatorWith(request: Mutator<any, any, any>["request"]) {
    const restore = vi.fn();
    const onSuccess = vi.fn();
    const mutator: Mutator<any, any, any> = {
      kind: "set_exposure_status",
      onMutate: vi.fn(() => ({ restore })),
      request,
      onSuccess,
    };
    return { mutator, restore, onSuccess };
  }

  it("stays in storage after persistence attaches, so a second reload keeps it too", () => {
    storeOps("op-a");
    const { mutator } = mutatorWith(downFor(60_000));
    void rehydrate(qc, new Map<OpKind, Mutator<any, any, any>>([["set_exposure_status", mutator]]));
    attachPersistence(qc.getMutationCache());
    expect(stored().map((o) => o.clientOpId)).toEqual(["op-a"]);
  });

  it("is pending in the mutation cache, where the outage banner can see it", () => {
    storeOps("op-a");
    const { mutator } = mutatorWith(downFor(60_000));
    void rehydrate(qc, new Map([["set_exposure_status" as OpKind, mutator]]));
    const pending = qc.getMutationCache().getAll().filter((m) => m.state.status === "pending");
    expect(pending).toHaveLength(1);
  });

  it("lands once the server is back, with the same client op id and no rollback", async () => {
    storeOps("op-a");
    const request = downFor(90_000);
    const { mutator, restore, onSuccess } = mutatorWith(request);
    const done = rehydrate(qc, new Map([["set_exposure_status" as OpKind, mutator]]));
    attachPersistence(qc.getMutationCache());

    await vi.advanceTimersByTimeAsync(60_000);
    expect(onSuccess).not.toHaveBeenCalled();
    expect(restore).not.toHaveBeenCalled();

    await vi.advanceTimersByTimeAsync(40_000);
    await expect(done).resolves.toEqual({ dropped: 0, replayed: 1, failed: 0, lost: 0 });
    expect(onSuccess).toHaveBeenCalledTimes(1);
    expect(restore).not.toHaveBeenCalled();
    const ids = request.mock.calls.map((c) => (c[0] as { clientOpId: string }).clientOpId);
    expect(new Set(ids)).toEqual(new Set(["op-a"]));
    expect(stored()).toEqual([]);
  });

  it("sends replayed edits one at a time, in the order they were made", async () => {
    storeOps("op-1", "op-2", "op-3");
    let release!: () => void;
    const held = new Promise<void>((r) => { release = r; });
    const request = vi.fn(async (p: { clientOpId: string }) => {
      if (p.clientOpId === "op-1") await held;
      return { ok: true };
    });
    const { mutator } = mutatorWith(request);
    const done = rehydrate(qc, new Map([["set_exposure_status" as OpKind, mutator]]));
    await vi.advanceTimersByTimeAsync(1_000);
    expect(request.mock.calls.map((c) => (c[0] as { clientOpId: string }).clientOpId)).toEqual(["op-1"]);
    release();
    await done;
    expect(request.mock.calls.map((c) => (c[0] as { clientOpId: string }).clientOpId)).toEqual(["op-1", "op-2", "op-3"]);
  });

  it("does not replay an op that is already pending in this page", () => {
    storeOps("op-a");
    const request = downFor(60_000);
    const { mutator } = mutatorWith(request);
    const mutators = new Map([["set_exposure_status" as OpKind, mutator]]);
    void rehydrate(qc, mutators);
    attachPersistence(qc.getMutationCache());
    void rehydrate(qc, mutators);
    expect(qc.getMutationCache().getAll().filter((m) => m.state.status === "pending")).toHaveLength(1);
    expect(mutator.onMutate).toHaveBeenCalledTimes(1);
  });

  it("skips a malformed stored entry and still replays the rest", async () => {
    sessionStorage.setItem(STORAGE_KEY, JSON.stringify([null, {
      schemaVersion: SCHEMA_VERSION, kind: "set_exposure_status", clientOpId: "op-b",
      payload: { kind: "set_exposure_status", clientOpId: "op-b", exposureId: 7, status: "rejected" },
    }]));
    const { mutator, onSuccess } = mutatorWith(vi.fn(async () => ({ ok: true })));
    const result = await rehydrate(qc, new Map([["set_exposure_status" as OpKind, mutator]]));
    expect(result).toEqual({ dropped: 0, replayed: 1, failed: 0, lost: 1 });
    expect(onSuccess).toHaveBeenCalledTimes(1);
  });

  it("counts an edit this build can no longer apply as lost, apart from one the server rejected", async () => {
    storeOps("op-stale", "op-rejected");
    const request = vi.fn(async () => {
      throw Object.assign(new Error("rejected"), { status: 422 });
    });
    const { mutator } = mutatorWith(request);
    mutator.onMutate = vi.fn((payload: { clientOpId: string }) => {
      if (payload.clientOpId === "op-stale") throw new Error("shape changed");
      return { restore: vi.fn() };
    });
    const result = await rehydrate(qc, new Map([["set_exposure_status" as OpKind, mutator]]));
    expect(result).toEqual({ dropped: 0, replayed: 0, failed: 1, lost: 1 });
    expect(request).toHaveBeenCalledTimes(1);
  });
});
