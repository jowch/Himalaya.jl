import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useSpeculativeSnap,
  useCreateSpeculative,
  queryKeys,
} from "../src/queries";
import { pendingDeferreds } from "../src/lib/queue/deferred";
import { makeFakeMutation } from "./queue/helpers";
import type { IndexEntry, SpeculativeSnap } from "../src/api";

const EXPOSURE_ID = 42;

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

function makeIndex(id: number): IndexEntry {
  return {
    id,
    exposure_id: EXPOSURE_ID,
    phase: "Pn3m",
    basis: 100,
    score: 0.9,
    r_squared: 0.99,
    lattice_d: 50,
    ngc: null,
    status: "candidate",
    kind: "speculative",
    inputs_hash: null,
    peaks: [],
    predicted_q: [],
  };
}


function makeSnap(): SpeculativeSnap[] {
  return [
    {
      ratio_position: 0,
      ratio_value: 1.0,
      target_q: 0.1,
      suggested_peak_id: 5,
      suggested_q: 0.1,
      suggested_residual: 0.0,
      is_anchor: true,
    },
  ];
}

describe("useSpeculativeSnap gating (M2.4)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
  });

  it("query enabled when no peak ops are pending", async () => {
    const { wrapper } = withClient();
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValueOnce(
      new Response(JSON.stringify(makeSnap()), {
        status: 200, headers: { "Content-Type": "application/json" },
      }),
    );
    const { result } = renderHook(
      () => useSpeculativeSnap(EXPOSURE_ID, "Pn3m", 5, 1.0),
      { wrapper },
    );
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(fetchSpy).toHaveBeenCalled();
    const url = (fetchSpy.mock.calls[0]?.[0] ?? "").toString();
    expect(url).toContain(`/api/exposures/${EXPOSURE_ID}/speculative-snap`);
  });

  it("query is disabled while a peak_added op is pending for the same exposure", () => {
    const { client, wrapper } = withClient();
    client.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_added", exposureId: EXPOSURE_ID, clientOpId: "op-1", payload: { q: 1.0 } },
    }));
    const fetchSpy = vi.spyOn(global, "fetch");
    const { result } = renderHook(
      () => useSpeculativeSnap(EXPOSURE_ID, "Pn3m", 5, 1.0),
      { wrapper },
    );
    expect(result.current.fetchStatus).toBe("idle");
    expect(result.current.isSuccess).toBe(false);
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it("query is disabled while a peak_excluded op is pending", () => {
    const { client, wrapper } = withClient();
    client.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_excluded", exposureId: EXPOSURE_ID, clientOpId: "op-x", payload: { peakId: 5 } },
    }));
    const fetchSpy = vi.spyOn(global, "fetch");
    const { result } = renderHook(
      () => useSpeculativeSnap(EXPOSURE_ID, "Pn3m", 5, 1.0),
      { wrapper },
    );
    expect(result.current.fetchStatus).toBe("idle");
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it("query is NOT disabled by a non-peak-affecting op (e.g. add_tag) on the same exposure", async () => {
    const { client, wrapper } = withClient();
    client.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "add_tag", exposureId: EXPOSURE_ID, clientOpId: "op-tag", payload: { key: "k", value: "v" } },
    }));
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValueOnce(
      new Response(JSON.stringify(makeSnap()), {
        status: 200, headers: { "Content-Type": "application/json" },
      }),
    );
    const { result } = renderHook(
      () => useSpeculativeSnap(EXPOSURE_ID, "Pn3m", 5, 1.0),
      { wrapper },
    );
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(fetchSpy).toHaveBeenCalled();
  });

  it("query stays enabled when a pending peak op is for a DIFFERENT exposure", async () => {
    const { client, wrapper } = withClient();
    client.getMutationCache().add(makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_added", exposureId: 99, clientOpId: "op-other", payload: { q: 1.0 } },
    }));
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValueOnce(
      new Response(JSON.stringify(makeSnap()), {
        status: 200, headers: { "Content-Type": "application/json" },
      }),
    );
    const { result } = renderHook(
      () => useSpeculativeSnap(EXPOSURE_ID, "Pn3m", 5, 1.0),
      { wrapper },
    );
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(fetchSpy).toHaveBeenCalled();
  });

  it("query is disabled when phase or anchorPeakId is undefined (independent of pending peak ops)", () => {
    const { wrapper } = withClient();
    const fetchSpy = vi.spyOn(global, "fetch");
    const { result } = renderHook(
      () => useSpeculativeSnap(EXPOSURE_ID, undefined, undefined, 1.0),
      { wrapper },
    );
    expect(result.current.fetchStatus).toBe("idle");
    expect(fetchSpy).not.toHaveBeenCalled();
  });

  it("query re-enables after the pending peak op transitions out of 'pending'", async () => {
    const { client, wrapper } = withClient();
    const mutation = makeFakeMutation({
      status: "pending",
      variables: { kind: "peak_added", exposureId: EXPOSURE_ID, clientOpId: "op-1", payload: { q: 1.0 } },
    });
    const mc = client.getMutationCache();
    mc.add(mutation);
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValueOnce(
      new Response(JSON.stringify(makeSnap()), {
        status: 200, headers: { "Content-Type": "application/json" },
      }),
    );
    const { result } = renderHook(
      () => useSpeculativeSnap(EXPOSURE_ID, "Pn3m", 5, 1.0),
      { wrapper },
    );
    expect(result.current.fetchStatus).toBe("idle");
    expect(fetchSpy).not.toHaveBeenCalled();
    // Flip the mutation status to success and notify mutation cache subscribers.
    act(() => {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      (mutation.state as any).status = "success";
      mc.notify({ type: "updated", mutation, action: { type: "success", data: undefined } } as never);
    });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(fetchSpy).toHaveBeenCalled();
  });
});

describe("useCreateSpeculative — createSpeculativeMutator (M2.4)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
  });

  it("calls api.createSpeculative with the mutator-built body shape", async () => {
    const { wrapper } = withClient();
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValueOnce(
      new Response(JSON.stringify(makeIndex(101)), {
        status: 200, headers: { "Content-Type": "application/json" },
      }),
    );
    const { result } = renderHook(() => useCreateSpeculative(EXPOSURE_ID), { wrapper });
    act(() => {
      result.current.mutate({
        phase: "Pn3m",
        anchor_peak_id: 5,
        anchor_ratio: 1.0,
        additional: [{ ratio_position: 1, peak_id: 6 }],
        active: true,
      });
    });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(fetchSpy).toHaveBeenCalled();
    const url = (fetchSpy.mock.calls[0]?.[0] ?? "").toString();
    expect(url).toBe(`/api/exposures/${EXPOSURE_ID}/speculative`);
    const init = fetchSpy.mock.calls[0]?.[1] as RequestInit | undefined;
    expect(init?.method).toBe("POST");
    const body = JSON.parse(String(init?.body ?? "{}"));
    expect(body).toEqual({
      phase: "Pn3m",
      anchor_peak_id: 5,
      anchor_ratio: 1.0,
      additional: [{ ratio_position: 1, peak_id: 6 }],
      active: true,
    });
  });

  it("threads X-Client-Op-Id header on the request", async () => {
    const { wrapper } = withClient();
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValueOnce(
      new Response(JSON.stringify(makeIndex(102)), {
        status: 200, headers: { "Content-Type": "application/json" },
      }),
    );
    const { result } = renderHook(() => useCreateSpeculative(EXPOSURE_ID), { wrapper });
    act(() => {
      result.current.mutate({
        phase: "Pn3m",
        anchor_peak_id: 5,
        anchor_ratio: 1.0,
        additional: [],
      });
    });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    const init = fetchSpy.mock.calls[0]?.[1] as RequestInit | undefined;
    const headers = new Headers(init?.headers);
    expect(headers.get("X-Client-Op-Id")).toBeTruthy();
    expect(headers.get("X-Client-Id")).toBeTruthy();
  });

  it("HTTP success splices the new index into the indices cache", async () => {
    const { client, wrapper } = withClient();
    const seedIndex = makeIndex(50);
    client.setQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID), [seedIndex]);
    mockOnce(200, makeIndex(101));
    const { result } = renderHook(() => useCreateSpeculative(EXPOSURE_ID), { wrapper });
    act(() => {
      result.current.mutate({
        phase: "Pn3m",
        anchor_peak_id: 5,
        anchor_ratio: 1.0,
        additional: [],
      });
    });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    const list = client.getQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID)) ?? [];
    expect(list.map((ix) => ix.id).sort((a, b) => a - b)).toEqual([50, 101]);
  });

  it("HTTP success is idempotent: replaying the response does not duplicate the index", async () => {
    const { client, wrapper } = withClient();
    const existing = makeIndex(101);
    client.setQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID), [existing]);
    mockOnce(200, makeIndex(101));
    const { result } = renderHook(() => useCreateSpeculative(EXPOSURE_ID), { wrapper });
    act(() => {
      result.current.mutate({
        phase: "Pn3m",
        anchor_peak_id: 5,
        anchor_ratio: 1.0,
        additional: [],
      });
    });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    const list = client.getQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID)) ?? [];
    expect(list).toHaveLength(1);
    expect(list[0]?.id).toBe(101);
  });

  it("HTTP error rolls back: indices cache restored to pre-mutate snapshot", async () => {
    const { client, wrapper } = withClient();
    const seedIndices = [makeIndex(50)];
    client.setQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID), seedIndices);
    mockOnce(400, { error: "bad anchor" });
    const { result } = renderHook(() => useCreateSpeculative(EXPOSURE_ID), { wrapper });
    act(() => {
      result.current.mutate({
        phase: "Pn3m",
        anchor_peak_id: 5,
        anchor_ratio: 1.0,
        additional: [],
      });
    });
    await waitFor(() => expect(result.current.error).toBeTruthy());
    const indices = client.getQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID));
    expect(indices).toEqual(seedIndices);
  });

  it("optimistic phase leaves caches unchanged (mutator is null-optimistic)", async () => {
    const { client, wrapper } = withClient();
    const seedIndices = [makeIndex(50)];
    client.setQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID), seedIndices);
    // Hold the request open so we can observe state during the pending phase.
    vi.spyOn(global, "fetch").mockImplementation(() => new Promise(() => {}));
    const { result } = renderHook(() => useCreateSpeculative(EXPOSURE_ID), { wrapper });
    act(() => {
      result.current.mutate({
        phase: "Pn3m",
        anchor_peak_id: 5,
        anchor_ratio: 1.0,
        additional: [],
      });
    });
    await waitFor(() => expect(result.current.isPending).toBe(true));
    const list = client.getQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID));
    expect(list).toEqual(seedIndices);
  });
});
