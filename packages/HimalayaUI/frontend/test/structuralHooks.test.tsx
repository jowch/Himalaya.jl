import { describe, it, expect } from "vitest";
import { renderHook } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import {
  useRenameSample, useMoveExposure, useMergeSamples, useSplitSample,
  useDismissGroupingFlag, useUpdateExperiment,
} from "../src/queries";

function wrapper({ children }: { children: ReactNode }) {
  const qc = new QueryClient();
  return <QueryClientProvider client={qc}>{children}</QueryClientProvider>;
}

describe("structural-edit + experiment hooks expose mutate()", () => {
  // All row-scoped hooks take only experimentId; entity ids (sampleId/exposureId)
  // go in the mutate() input (one hook instance per experiment, not per row).
  it("useRenameSample", () => { expect(typeof renderHook(() => useRenameSample(7), { wrapper }).result.current.mutate).toBe("function"); });
  it("useMoveExposure", () => { expect(typeof renderHook(() => useMoveExposure(7), { wrapper }).result.current.mutate).toBe("function"); });
  it("useMergeSamples", () => { expect(typeof renderHook(() => useMergeSamples(7), { wrapper }).result.current.mutate).toBe("function"); });
  it("useSplitSample", () => { expect(typeof renderHook(() => useSplitSample(7), { wrapper }).result.current.mutate).toBe("function"); });
  it("useDismissGroupingFlag", () => { expect(typeof renderHook(() => useDismissGroupingFlag(7), { wrapper }).result.current.mutate).toBe("function"); });
  it("useUpdateExperiment", () => { expect(typeof renderHook(() => useUpdateExperiment(7), { wrapper }).result.current.mutate).toBe("function"); });
});
