import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { useScreenedProgress, queryKeys } from "../src/queries";
import type { CorpusSample, Exposure } from "../src/api";

// SA-SCREENED integration pin: the contact-sheet hero metric ("N / M samples
// screened") must advance once every exposure of a sample is ACCEPTED — the
// Keep verb, not just Drop, completes screening. The per-sample derivation is
// covered in screened.test.ts; this pins the aggregate the header renders.

const SAMPLE: CorpusSample = {
  id: 1, experiment_id: 1, name: "s1", display_name: null, notes: "",
  tags: [], q_units: "A-1",
} as CorpusSample;

function exposure(id: number, status: "accepted" | "rejected" | null): Exposure {
  return {
    id, sample_id: 1, filename: `f${id}.dat`, kind: "file", selected: false,
    status, image_path: null, image_version: "", tags: [], sources: [],
    trace_hash: null, analysis_inputs_hash: null,
  };
}

function setup(rows: Exposure[]) {
  const client = makeClient();
  client.setQueryData(queryKeys.exposures(1), rows);
  // The seeded rows are stale (staleTime 0) → useQueries refetches; answer
  // with the same rows so the assertion is race-free either way.
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify(rows), {
      status: 200, headers: { "Content-Type": "application/json" },
    }),
  );
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return renderHook(() => useScreenedProgress([SAMPLE]), { wrapper });
}

describe("useScreenedProgress — the Keep verb completes screening", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("counts a sample as screened when every exposure is accepted", () => {
    const { result } = setup([exposure(10, "accepted"), exposure(11, "accepted")]);
    expect(result.current).toEqual({ screened: 1, total: 1 });
  });

  it("a mix of accepted and rejected still screens (every status non-null)", () => {
    const { result } = setup([exposure(10, "accepted"), exposure(11, "rejected")]);
    expect(result.current).toEqual({ screened: 1, total: 1 });
  });

  it("one unscreened (null) exposure keeps the sample out of the count", () => {
    const { result } = setup([exposure(10, "accepted"), exposure(11, null)]);
    expect(result.current).toEqual({ screened: 0, total: 1 });
  });
});
