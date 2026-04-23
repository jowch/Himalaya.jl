import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useUpdateSample, useAddSampleTag, useRemoveSampleTag, queryKeys,
} from "../src/queries";

const EXPERIMENT_ID = 1;

function withClient() {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, wrapper };
}

function mockOnce(status: number, body: unknown): void {
  vi.spyOn(global, "fetch").mockResolvedValueOnce(
    new Response(status === 204 ? null : JSON.stringify(body), {
      status, headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("queries — sample mutations", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("useUpdateSample invalidates samples on success", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(200, { id: 10, experiment_id: 1, label: "A1", name: "s1",
                    notes: "updated", tags: [] });
    const { result } = renderHook(
      () => useUpdateSample(EXPERIMENT_ID, 10), { wrapper },
    );
    await act(async () => { await result.current.mutateAsync({ notes: "updated" }); });
    expect(invalidate).toHaveBeenCalledWith({
      queryKey: queryKeys.samples(EXPERIMENT_ID),
    });
  });

  it("useAddSampleTag invalidates samples on success", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(201, { id: 5, sample_id: 10, key: "lipid", value: "DOPC",
                    source: "manual" });
    const { result } = renderHook(
      () => useAddSampleTag(EXPERIMENT_ID, 10), { wrapper },
    );
    await act(async () => {
      await result.current.mutateAsync({ key: "lipid", value: "DOPC" });
    });
    expect(invalidate).toHaveBeenCalledWith({
      queryKey: queryKeys.samples(EXPERIMENT_ID),
    });
  });

  it("useRemoveSampleTag invalidates samples on success", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(204, null);
    const { result } = renderHook(
      () => useRemoveSampleTag(EXPERIMENT_ID, 10), { wrapper },
    );
    await act(async () => { await result.current.mutateAsync(5); });
    expect(invalidate).toHaveBeenCalledWith({
      queryKey: queryKeys.samples(EXPERIMENT_ID),
    });
  });
});
