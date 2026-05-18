import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { useCorpusSamples, queryKeys } from "../src/queries";

function withClient() {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, wrapper };
}

function mockOnce(status: number, body: unknown): void {
  vi.spyOn(global, "fetch").mockResolvedValueOnce(
    new Response(JSON.stringify(body), {
      status, headers: { "Content-Type": "application/json" },
    }),
  );
}

const CORPUS_ROWS = [
  { id: 1, experiment_id: 1, name: "s1", display_name: null, notes: "",
    tags: [], q_units: "A-1" },
  { id: 2, experiment_id: 2, name: "s2", display_name: null, notes: "",
    tags: [{ id: 9, key: "lipid", value: "DOPC", source: "manual" }],
    q_units: "nm-1" },
];

describe("queries — useCorpusSamples", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("fetches GET /api/samples and returns the corpus list with q_units", async () => {
    const { wrapper } = withClient();
    mockOnce(200, CORPUS_ROWS);
    const { result } = renderHook(() => useCorpusSamples(), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));

    expect(global.fetch).toHaveBeenCalledWith(
      "/api/samples", expect.objectContaining({ method: "GET" }),
    );
    expect(result.current.data).toHaveLength(2);
    expect(result.current.data?.[0].q_units).toBe("A-1");
    expect(result.current.data?.[1].q_units).toBe("nm-1");
    expect(result.current.data?.[1].tags[0].key).toBe("lipid");
  });

  it("uses the ['corpus','samples'] key, distinct from per-experiment samples", () => {
    expect(queryKeys.corpusSamples).toEqual(["corpus", "samples"]);
    // The per-experiment key is experiment-scoped for every id — the corpus
    // key must never deep-equal it, so cache entries cannot collide.
    for (const id of [0, 1, 42]) {
      expect(queryKeys.corpusSamples).not.toEqual(queryKeys.samples(id));
    }
  });
});
