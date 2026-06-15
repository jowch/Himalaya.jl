/**
 * useSetExposureStatusBatch — cross-sample cull dispatch hook.
 *
 * The key invariant: useQueueMutation merges payload as
 *   { kind, clientOpId, ...scope, ...input }
 * so a sampleId passed in the per-call INPUT overrides the placeholder
 * scope.sampleId (0) and lands in the mutator as the correct sample.
 *
 * We mock useQueueMutation so the inner mutate is a spy, then assert that
 * the batch wrapper forwards { sampleId, exposureId, status } untouched.
 * useAppState is mocked to provide a stable username without Zustand setup.
 */
import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { useSetExposureStatusBatch } from "../src/queries";

// ---------------------------------------------------------------------------
// Module mocks
// ---------------------------------------------------------------------------

const innerMutate = vi.fn();

vi.mock("../src/lib/queue/useQueueMutation", () => ({
  useQueueMutation: vi.fn(() => ({
    mutate: innerMutate,
    isPending: false,
    isSuccess: false,
    data: undefined,
    error: null,
    reset: vi.fn(),
  })),
}));

vi.mock("../src/state", () => ({
  useAppState: vi.fn((selector: (s: { username: string }) => unknown) =>
    selector({ username: "testuser" }),
  ),
}));

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function withClient() {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { wrapper };
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

describe("useSetExposureStatusBatch", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    // Re-initialise innerMutate so each test gets a clean spy.
    innerMutate.mockReset();
  });

  it("forwards sampleId from the per-call input to the inner mutate", () => {
    const { wrapper } = withClient();
    const { result } = renderHook(() => useSetExposureStatusBatch(), { wrapper });

    act(() => {
      result.current.mutate({ sampleId: 2, exposureId: 99, status: "rejected" });
    });

    expect(innerMutate).toHaveBeenCalledOnce();
    const calledWith = innerMutate.mock.calls[0][0] as Record<string, unknown>;
    expect(calledWith).toMatchObject({ sampleId: 2, exposureId: 99, status: "rejected" });
  });

  it("accepts status: null (restore) with the correct sampleId", () => {
    const { wrapper } = withClient();
    const { result } = renderHook(() => useSetExposureStatusBatch(), { wrapper });

    act(() => {
      result.current.mutate({ sampleId: 5, exposureId: 7, status: null });
    });

    expect(innerMutate).toHaveBeenCalledOnce();
    const calledWith = innerMutate.mock.calls[0][0] as Record<string, unknown>;
    expect(calledWith).toMatchObject({ sampleId: 5, exposureId: 7, status: null });
  });

  it("accepts status: accepted with the correct sampleId", () => {
    const { wrapper } = withClient();
    const { result } = renderHook(() => useSetExposureStatusBatch(), { wrapper });

    act(() => {
      result.current.mutate({ sampleId: 99, exposureId: 12, status: "accepted" });
    });

    expect(innerMutate).toHaveBeenCalledOnce();
    const calledWith = innerMutate.mock.calls[0][0] as Record<string, unknown>;
    expect(calledWith).toMatchObject({ sampleId: 99, exposureId: 12, status: "accepted" });
  });

  it("each call is independent (different sampleIds dispatch correctly)", () => {
    const { wrapper } = withClient();
    const { result } = renderHook(() => useSetExposureStatusBatch(), { wrapper });

    act(() => {
      result.current.mutate({ sampleId: 1, exposureId: 10, status: "rejected" });
      result.current.mutate({ sampleId: 2, exposureId: 20, status: "accepted" });
      result.current.mutate({ sampleId: 3, exposureId: 30, status: null });
    });

    expect(innerMutate).toHaveBeenCalledTimes(3);
    expect(innerMutate.mock.calls[0][0]).toMatchObject({ sampleId: 1, exposureId: 10, status: "rejected" });
    expect(innerMutate.mock.calls[1][0]).toMatchObject({ sampleId: 2, exposureId: 20, status: "accepted" });
    expect(innerMutate.mock.calls[2][0]).toMatchObject({ sampleId: 3, exposureId: 30, status: null });
  });

  it("exposes isPending, isSuccess, error, reset from the inner hook result", () => {
    const { wrapper } = withClient();
    const { result } = renderHook(() => useSetExposureStatusBatch(), { wrapper });

    // Defaults from the mock.
    expect(result.current.isPending).toBe(false);
    expect(result.current.isSuccess).toBe(false);
    expect(result.current.error).toBeNull();
    expect(typeof result.current.reset).toBe("function");
  });
});
