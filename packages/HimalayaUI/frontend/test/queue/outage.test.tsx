import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { renderHook, act } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { useQueueMutation } from "../../src/lib/queue/useQueueMutation";
import type { Mutator } from "../../src/lib/queue/types";
import { pendingDeferreds } from "../../src/lib/queue/deferred";
import * as toast from "../../src/lib/toast";
import { NetworkError } from "../../src/api";

function wrap(qc: QueryClient) {
  return ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={qc}>{children}</QueryClientProvider>
  );
}

function serverDownFor(ms: number, response: unknown = { ok: true }) {
  const upAt = Date.now() + ms;
  return vi.fn(async () => {
    if (Date.now() < upAt) throw new NetworkError("PUT", "/api/x");
    return response;
  });
}

describe("an edit made while the server is down", () => {
  let qc: QueryClient;
  let toastSpy: ReturnType<typeof vi.spyOn>;

  beforeEach(() => {
    vi.useFakeTimers();
    qc = new QueryClient();
    pendingDeferreds.clear();
    toastSpy = vi.spyOn(toast, "showToast").mockImplementation(() => {});
  });
  afterEach(() => {
    toastSpy.mockRestore();
    vi.useRealTimers();
  });

  it("is applied once the server comes back after a restart longer than a minute", async () => {
    const restore = vi.fn();
    const onSuccess = vi.fn();
    const mutator: Mutator<{ status: string }, Record<string, unknown>, unknown> = {
      kind: "set_exposure_status",
      onMutate: vi.fn(() => ({ restore })),
      request: serverDownFor(90_000),
      onSuccess,
    };
    const { result } = renderHook(() => useQueueMutation(mutator, { exposureId: 7 }), {
      wrapper: wrap(qc),
    });
    const onError = vi.fn();
    act(() => result.current.mutate({ status: "rejected" }, { onError }));

    await act(async () => { await vi.advanceTimersByTimeAsync(60_000); });
    expect(result.current.isPending).toBe(true);
    expect(restore).not.toHaveBeenCalled();

    await act(async () => { await vi.advanceTimersByTimeAsync(40_000); });
    expect(onSuccess).toHaveBeenCalledTimes(1);
    expect(restore).not.toHaveBeenCalled();
    expect(onError).not.toHaveBeenCalled();
    const ids = (mutator.request as ReturnType<typeof vi.fn>).mock.calls.map(
      (c) => (c[0] as { clientOpId: string }).clientOpId,
    );
    expect(new Set(ids).size).toBe(1);
  });

  it("retries at most every 5 s while the server is down, so it lands soon after the server is back", async () => {
    const onSuccess = vi.fn();
    const mutator: Mutator<{ status: string }, Record<string, unknown>, unknown> = {
      kind: "set_exposure_status",
      onMutate: vi.fn(() => ({ restore: vi.fn() })),
      request: serverDownFor(120_000),
      onSuccess,
    };
    const { result } = renderHook(() => useQueueMutation(mutator, {}), { wrapper: wrap(qc) });
    act(() => result.current.mutate({ status: "rejected" }));
    await act(async () => { await vi.advanceTimersByTimeAsync(120_000); });
    expect(onSuccess).not.toHaveBeenCalled();
    await act(async () => { await vi.advanceTimersByTimeAsync(5_000); });
    expect(onSuccess).toHaveBeenCalledTimes(1);
  });

  it("says so when a server error outlasts the retries, instead of undoing it silently", async () => {
    const restore = vi.fn();
    const mutator: Mutator<{ status: string }, Record<string, unknown>, unknown> = {
      kind: "set_exposure_status",
      onMutate: vi.fn(() => ({ restore })),
      request: vi.fn(async () => {
        throw Object.assign(new Error("500 Internal Server Error"), { status: 500 });
      }),
      onSuccess: vi.fn(),
    };
    const { result } = renderHook(() => useQueueMutation(mutator, {}), { wrapper: wrap(qc) });
    act(() => result.current.mutate({ status: "rejected" }));
    await act(async () => { await vi.advanceTimersByTimeAsync(120_000); });
    expect(restore).toHaveBeenCalledTimes(1);
    expect(toastSpy).toHaveBeenCalledWith(
      expect.stringMatching(/^Couldn't update exposure status/),
      "error",
    );
  });

  it("delivers a sample's queued edits in the order they were made", async () => {
    const upAt = Date.now() + 60_000;
    const applied: string[] = [];
    const mutator: Mutator<{ status: string }, Record<string, unknown>, unknown> = {
      kind: "set_exposure_status",
      onMutate: () => ({ restore: () => {} }),
      request: async (p) => {
        if (Date.now() < upAt) throw new NetworkError("PUT", "/api/x");
        applied.push((p as unknown as { status: string }).status);
        return {};
      },
      onSuccess: () => {},
    };
    const { result } = renderHook(() => useQueueMutation(mutator, { sampleId: 1 }), { wrapper: wrap(qc) });
    act(() => result.current.mutate({ status: "rejected" }));
    await act(async () => { await vi.advanceTimersByTimeAsync(3_000); });
    act(() => result.current.mutate({ status: "accepted" }));
    await act(async () => { await vi.advanceTimersByTimeAsync(80_000); });
    expect(applied).toEqual(["rejected", "accepted"]);
  });

  it("does not hold an edit on one sample behind a stuck edit on another", async () => {
    const landed: number[] = [];
    const stuck: Mutator<{ status: string }, Record<string, unknown>, unknown> = {
      kind: "set_exposure_status",
      onMutate: () => ({ restore: () => {} }),
      request: async () => { throw new NetworkError("PUT", "/api/a"); },
      onSuccess: () => {},
    };
    const fine: Mutator<{ status: string }, Record<string, unknown>, unknown> = {
      kind: "set_exposure_status",
      onMutate: () => ({ restore: () => {} }),
      request: async () => { landed.push(2); return {}; },
      onSuccess: () => {},
    };
    const a = renderHook(() => useQueueMutation(stuck, { sampleId: 1 }), { wrapper: wrap(qc) });
    const b = renderHook(() => useQueueMutation(fine, { sampleId: 2 }), { wrapper: wrap(qc) });
    act(() => a.result.current.mutate({ status: "rejected" }));
    act(() => b.result.current.mutate({ status: "rejected" }));
    await act(async () => { await vi.advanceTimersByTimeAsync(1_000); });
    expect(landed).toEqual([2]);
    expect(a.result.current.isPending).toBe(true);
  });

  it("gives a client-side throw the bounded budget and no server-error toast", async () => {
    const restore = vi.fn();
    const mutator: Mutator<{ status: string }, Record<string, unknown>, unknown> = {
      kind: "set_exposure_status",
      onMutate: vi.fn(() => ({ restore })),
      request: vi.fn(async () => { throw new TypeError("Cannot read properties of undefined"); }),
      onSuccess: vi.fn(),
    };
    const { result } = renderHook(() => useQueueMutation(mutator, {}), { wrapper: wrap(qc) });
    act(() => result.current.mutate({ status: "rejected" }));
    await act(async () => { await vi.advanceTimersByTimeAsync(120_000); });
    expect(restore).toHaveBeenCalledTimes(1);
    expect(toastSpy).not.toHaveBeenCalled();
  });

  it("leaves the error message to a caller that handles its own failure", async () => {
    const mutator: Mutator<{ status: string }, Record<string, unknown>, unknown> = {
      kind: "set_exposure_status",
      onMutate: () => ({ restore: () => {} }),
      request: async () => { throw Object.assign(new Error("500"), { status: 500 }); },
      onSuccess: () => {},
    };
    const onError = vi.fn();
    const { result } = renderHook(() => useQueueMutation(mutator, {}), { wrapper: wrap(qc) });
    act(() => result.current.mutate({ status: "rejected" }, { onError }));
    await act(async () => { await vi.advanceTimersByTimeAsync(120_000); });
    expect(onError).toHaveBeenCalledTimes(1);
    expect(toastSpy).not.toHaveBeenCalled();
  });
});
