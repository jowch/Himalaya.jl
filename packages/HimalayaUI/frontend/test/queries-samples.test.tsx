import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useUpdateSample, useAddSampleTag, useRemoveSampleTag, queryKeys,
} from "../src/queries";
import { pendingDeferreds } from "../src/lib/queue/deferred";

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

describe("queries — sample mutations (queue-driven)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
  });

  it("useUpdateSample writes the server response into the samples cache on success", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData(queryKeys.samples(EXPERIMENT_ID), [
      { id: 10, experiment_id: 1, label: "A1", name: "s1", notes: "", tags: [] },
    ]);
    mockOnce(200, { id: 10, experiment_id: 1, label: "A1", name: "s1",
                    notes: "updated", tags: [] });
    const { result } = renderHook(
      () => useUpdateSample(EXPERIMENT_ID, 10), { wrapper },
    );
    act(() => { result.current.mutate({ notes: "updated" }); });
    await waitFor(() => {
      const list = client.getQueryData(queryKeys.samples(EXPERIMENT_ID)) as Array<{ notes: string }>;
      expect(list[0].notes).toBe("updated");
    });
  });

  it("useAddSampleTag appends an optimistic placeholder, replaced by server tag", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData(queryKeys.samples(EXPERIMENT_ID), [
      { id: 10, experiment_id: 1, label: "A1", name: "s1", notes: "", tags: [] },
    ]);
    mockOnce(201, { id: 5, sample_id: 10, key: "lipid", value: "DOPC",
                    source: "manual" });
    const { result } = renderHook(
      () => useAddSampleTag(EXPERIMENT_ID, 10), { wrapper },
    );
    act(() => { result.current.mutate({ key: "lipid", value: "DOPC" }); });
    await waitFor(() => {
      const list = client.getQueryData(queryKeys.samples(EXPERIMENT_ID)) as Array<{ tags: Array<{ id: number; key: string }> }>;
      expect(list[0].tags).toHaveLength(1);
      expect(list[0].tags[0].id).toBe(5);
      expect(list[0].tags[0].key).toBe("lipid");
    });
  });

  it("useRemoveSampleTag removes the tag optimistically", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData(queryKeys.samples(EXPERIMENT_ID), [
      { id: 10, experiment_id: 1, label: "A1", name: "s1", notes: "",
        tags: [{ id: 5, key: "lipid", value: "DOPC", source: "manual" }] },
    ]);
    mockOnce(204, null);
    const { result } = renderHook(
      () => useRemoveSampleTag(EXPERIMENT_ID, 10), { wrapper },
    );
    act(() => { result.current.mutate(5); });
    await waitFor(() => {
      const list = client.getQueryData(queryKeys.samples(EXPERIMENT_ID)) as Array<{ tags: unknown[] }>;
      expect(list[0].tags).toHaveLength(0);
    });
  });
});
