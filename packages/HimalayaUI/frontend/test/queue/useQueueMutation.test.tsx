import { describe, it, expect, beforeEach, vi } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { useQueueMutation } from "../../src/lib/queue/useQueueMutation";
import { handleRemoteEvent } from "../../src/lib/queue/replayCoordinator";
import type { Mutator, OpPayload } from "../../src/lib/queue/types";
import { pendingDeferreds } from "../../src/lib/queue/deferred";

function withQueryClient(qc: QueryClient) {
  return ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={qc}>{children}</QueryClientProvider>
  );
}

function makeMutator<TInput, TResponse>(
  partial: Partial<Mutator<OpPayload<TInput>, TResponse>> = {},
): Mutator<OpPayload<TInput>, TResponse> {
  return {
    kind: "peak_added",
    onMutate: vi.fn(() => ({ restore: vi.fn() })),
    request: vi.fn().mockResolvedValue({} as TResponse),
    onSuccess: vi.fn(),
    ...partial,
  };
}

describe("useQueueMutation", () => {
  let qc: QueryClient;
  beforeEach(() => {
    qc = new QueryClient({
      defaultOptions: { mutations: { retry: false } },
    });
    pendingDeferreds.clear();
  });

  it("mints a fresh clientOpId per mutate() call", async () => {
    const mutator = makeMutator();
    const { result } = renderHook(
      () => useQueueMutation(mutator, { exposureId: 42 }),
      { wrapper: withQueryClient(qc) },
    );

    act(() => result.current.mutate({ q: 1.0 }));
    act(() => result.current.mutate({ q: 2.0 }));

    await waitFor(() => expect((mutator.onMutate as any).mock.calls.length).toBe(2));
    const op1 = (mutator.onMutate as any).mock.calls[0][0] as OpPayload<unknown>;
    const op2 = (mutator.onMutate as any).mock.calls[1][0] as OpPayload<unknown>;
    expect(op1.clientOpId).not.toBe(op2.clientOpId);
    expect(op1.clientOpId).toMatch(/^[0-9a-f-]{36}$/);
  });

  it("merges scope + input into the payload passed to onMutate / request / onSuccess", async () => {
    const mutator = makeMutator<{ q: number }, { ok: true }>();
    (mutator.request as any).mockResolvedValue({ ok: true });
    const { result } = renderHook(
      () => useQueueMutation(mutator, { exposureId: 42, username: "alice" }),
      { wrapper: withQueryClient(qc) },
    );

    act(() => result.current.mutate({ q: 1.5 }));

    await waitFor(() => expect((mutator.onSuccess as any).mock.calls.length).toBe(1));
    const passed = (mutator.onMutate as any).mock.calls[0][0];
    expect(passed).toMatchObject({
      kind: "peak_added",
      exposureId: 42,
      username: "alice",
      q: 1.5,
    });
    expect(passed.clientOpId).toBeTruthy();
  });

  it("HTTP success → mutator.onSuccess called with the response", async () => {
    const mutator = makeMutator<{ q: number }, { id: number }>();
    (mutator.request as any).mockResolvedValue({ id: 99 });
    const { result } = renderHook(
      () => useQueueMutation(mutator, {}),
      { wrapper: withQueryClient(qc) },
    );

    act(() => result.current.mutate({ q: 3.0 }));
    await waitFor(() => expect((mutator.onSuccess as any).mock.calls.length).toBe(1));

    const [payload, response] = (mutator.onSuccess as any).mock.calls[0];
    expect(payload.q).toBe(3.0);
    expect(response).toEqual({ id: 99 });
  });

  it("HTTP 4xx error → context.restore called (rollback)", async () => {
    const restore = vi.fn();
    const mutator = makeMutator();
    (mutator.onMutate as any).mockReturnValue({ restore });
    (mutator.request as any).mockRejectedValue(
      Object.assign(new Error("400 bad"), { status: 400 }),
    );
    const { result } = renderHook(
      () => useQueueMutation(mutator, {}),
      { wrapper: withQueryClient(qc) },
    );

    act(() => result.current.mutate({ q: 1.0 } as any));
    await waitFor(() => expect(result.current.error).toBeTruthy());
    expect(restore).toHaveBeenCalledTimes(1);
  });

  it("validation error does NOT retry", async () => {
    const mutator = makeMutator();
    (mutator.request as any).mockRejectedValue(
      Object.assign(new Error("422"), { status: 422 }),
    );
    // Override the QueryClient to not block retry by default; let useQueueMutation's retry fn decide.
    qc = new QueryClient();
    const { result } = renderHook(
      () => useQueueMutation(mutator, {}),
      { wrapper: withQueryClient(qc) },
    );

    act(() => result.current.mutate({ q: 1.0 } as any));
    await waitFor(() => expect(result.current.error).toBeTruthy());
    // Only one call attempted (no retries).
    expect((mutator.request as any).mock.calls.length).toBe(1);
  });

  it("SSE event with matching client_op_id resolves before HTTP returns", async () => {
    let httpResolve: (v: any) => void = () => {};
    const slowHttp = new Promise((res) => { httpResolve = res; });
    const mutator = makeMutator();
    (mutator.request as any).mockReturnValue(slowHttp);

    const { result } = renderHook(
      () => useQueueMutation(mutator, {}),
      { wrapper: withQueryClient(qc) },
    );

    act(() => result.current.mutate({ q: 2.5 } as any));

    // Wait for the deferred to register.
    await waitFor(() => expect(pendingDeferreds.size).toBe(1));
    const opId = Array.from(pendingDeferreds.keys())[0];

    // Inject SSE before HTTP resolves.
    handleRemoteEvent({
      id: 7,
      kind: "peak_added",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: opId,
      payload: { q: 2.5 },
    }, qc, qc.getMutationCache());

    await waitFor(() => expect((mutator.onSuccess as any).mock.calls.length).toBe(1));
    expect((mutator.onSuccess as any).mock.calls[0][1]).toMatchObject({ event_id: 7 });

    // HTTP eventually resolves but should be a no-op (deferred already cleared).
    httpResolve({ id: 99 });
  });

  it("isPending reflects MutationCache state (true while in-flight, false after settle)", async () => {
    let httpResolve: (v: any) => void = () => {};
    const mutator = makeMutator();
    (mutator.request as any).mockImplementation(() => new Promise((res) => { httpResolve = res; }));
    const { result } = renderHook(
      () => useQueueMutation(mutator, {}),
      { wrapper: withQueryClient(qc) },
    );

    expect(result.current.isPending).toBe(false);
    act(() => result.current.mutate({ q: 1.0 } as any));
    await waitFor(() => expect(result.current.isPending).toBe(true));
    act(() => httpResolve({ id: 1 }));
    await waitFor(() => expect(result.current.isPending).toBe(false));
  });

  it("reset() clears error state", async () => {
    const mutator = makeMutator();
    (mutator.request as any).mockRejectedValue(
      Object.assign(new Error("400"), { status: 400 }),
    );
    const { result } = renderHook(
      () => useQueueMutation(mutator, {}),
      { wrapper: withQueryClient(qc) },
    );
    act(() => result.current.mutate({ q: 1.0 } as any));
    await waitFor(() => expect(result.current.error).toBeTruthy());
    act(() => result.current.reset());
    await waitFor(() => expect(result.current.error).toBeNull());
  });

  it("AbortSignal threads to mutator.request", async () => {
    let receivedSignal: AbortSignal | undefined;
    const mutator = makeMutator();
    (mutator.request as any).mockImplementation((_p: any, signal: AbortSignal) => {
      receivedSignal = signal;
      return new Promise(() => { /* never resolve */ });
    });
    const { result } = renderHook(
      () => useQueueMutation(mutator, {}),
      { wrapper: withQueryClient(qc) },
    );
    act(() => result.current.mutate({ q: 1.0 } as any));
    await waitFor(() => expect(receivedSignal).toBeDefined());
    expect(receivedSignal!.aborted).toBe(false);
  });
});
