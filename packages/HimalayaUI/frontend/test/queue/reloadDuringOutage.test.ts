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
  return vi.fn(async () => {
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
    await expect(done).resolves.toEqual({ dropped: 0, replayed: 1, failed: 0 });
    expect(onSuccess).toHaveBeenCalledTimes(1);
    expect(restore).not.toHaveBeenCalled();
    const ids = request.mock.calls.map((c) => (c[0] as { clientOpId: string }).clientOpId);
    expect(new Set(ids)).toEqual(new Set(["op-a"]));
    expect(stored()).toEqual([]);
  });

  it("lands replayed edits in the order they were made", async () => {
    storeOps("op-1", "op-2", "op-3");
    const landed: string[] = [];
    const upAt = Date.now() + 30_000;
    const request = vi.fn(async (p: { clientOpId: string }) => {
      if (Date.now() < upAt) throw new NetworkError("PUT", "/api/x");
      landed.push(p.clientOpId);
      return { ok: true };
    });
    const { mutator } = mutatorWith(request);
    const done = rehydrate(qc, new Map([["set_exposure_status" as OpKind, mutator]]));
    await vi.advanceTimersByTimeAsync(60_000);
    await done;
    expect(landed).toEqual(["op-1", "op-2", "op-3"]);
  });
});
