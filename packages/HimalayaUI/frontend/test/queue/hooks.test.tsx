import { describe, it, expect, beforeEach } from "vitest";
import { renderHook } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import {
  useExposureHasPendingPeakOps,
  useQueueOpStatus,
} from "../../src/lib/queue/hooks";
import { makeFakeMutation } from "./helpers";

function withQueryClient(qc: QueryClient) {
  return ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={qc}>{children}</QueryClientProvider>
  );
}

describe("useExposureHasPendingPeakOps", () => {
  let qc: QueryClient;
  beforeEach(() => {
    qc = new QueryClient();
  });

  it("returns false when no peak ops are pending for the exposure", () => {
    const { result } = renderHook(() => useExposureHasPendingPeakOps(42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe(false);
  });

  it("returns true while a peak_added op is pending for the same exposure", () => {
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_added", exposureId: 42, clientOpId: "op-1", payload: { q: 1.0 } },
    }));
    const { result } = renderHook(() => useExposureHasPendingPeakOps(42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe(true);
  });

  it("returns true while a reanalyze_exposure op is pending (peak-affecting)", () => {
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "reanalyze_exposure", exposureId: 42, clientOpId: "op-r", payload: {} },
    }));
    const { result } = renderHook(() => useExposureHasPendingPeakOps(42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe(true);
  });

  it("returns false for a different exposure id", () => {
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_added", exposureId: 99, clientOpId: "op-x", payload: { q: 1.0 } },
    }));
    const { result } = renderHook(() => useExposureHasPendingPeakOps(42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe(false);
  });

  it("returns false for a non-peak op kind on the same exposure", () => {
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "add_tag", exposureId: 42, clientOpId: "op-tag", payload: { key: "k", value: "v" } },
    }));
    const { result } = renderHook(() => useExposureHasPendingPeakOps(42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe(false);
  });

  it("returns false when exposureId is undefined", () => {
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_added", exposureId: 42, clientOpId: "op-1", payload: { q: 1.0 } },
    }));
    const { result } = renderHook(() => useExposureHasPendingPeakOps(undefined), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe(false);
  });
});

describe("useQueueOpStatus", () => {
  let qc: QueryClient;
  beforeEach(() => {
    qc = new QueryClient();
  });

  it("returns 'pending' when a matching kind is in flight (scoped by exposureId)", () => {
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "reanalyze_exposure", exposureId: 42, clientOpId: "op-r", payload: {} },
    }));
    const { result } = renderHook(() => useQueueOpStatus("reanalyze_exposure", 42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe("pending");
  });

  it("returns 'idle' for a non-matching kind", () => {
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_added", exposureId: 42, clientOpId: "op-1", payload: { q: 1.0 } },
    }));
    const { result } = renderHook(() => useQueueOpStatus("reanalyze_exposure", 42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe("idle");
  });

  it("ignores exposureId scope when not provided", () => {
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "reanalyze_exposure", exposureId: 999, clientOpId: "op-r", payload: {} },
    }));
    const { result } = renderHook(() => useQueueOpStatus("reanalyze_exposure"), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe("pending");
  });

  it("returns 'idle' for kind+exposureId mismatch", () => {
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "reanalyze_exposure", exposureId: 99, clientOpId: "op-r", payload: {} },
    }));
    const { result } = renderHook(() => useQueueOpStatus("reanalyze_exposure", 42), {
      wrapper: withQueryClient(qc),
    });
    expect(result.current).toBe("idle");
  });
});
