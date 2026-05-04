import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useAddIndexToGroup,
  useRemoveIndexFromGroup,
  useDeleteIndex,
  queryKeys,
} from "../src/queries";
import { pendingDeferreds } from "../src/lib/queue/deferred";
import type { GroupEntry, IndexEntry } from "../src/api";

const EXPOSURE_ID = 42;
const GROUP_ID = 7;

function makeGroup(members: number[] = []): GroupEntry {
  return {
    id: GROUP_ID,
    exposure_id: EXPOSURE_ID,
    kind: "custom",
    active: true,
    members,
  };
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

function mockNever(): void {
  vi.spyOn(global, "fetch").mockImplementation(() => new Promise(() => {}));
}

describe("queries — index/group mutations (queue-driven, M2.3)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
  });

  // ---------------------- useAddIndexToGroup ----------------------

  it("useAddIndexToGroup appends indexId to the group's members optimistically", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID), [makeGroup([])]);
    mockNever();
    const { result } = renderHook(() => useAddIndexToGroup(EXPOSURE_ID, GROUP_ID), { wrapper });
    act(() => { result.current.mutate(101); });
    await waitFor(() => {
      const groups = client.getQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID)) ?? [];
      expect(groups[0].members).toEqual([101]);
    });
  });

  it("useAddIndexToGroup replaces the optimistic group with the server-returned GroupEntry", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID), [makeGroup([])]);
    mockOnce(200, { ...makeGroup([101, 102]) });
    const { result } = renderHook(() => useAddIndexToGroup(EXPOSURE_ID, GROUP_ID), { wrapper });
    act(() => { result.current.mutate(101); });
    await waitFor(() => {
      const groups = client.getQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID)) ?? [];
      expect(groups[0].members).toEqual([101, 102]);
    });
  });

  it("useAddIndexToGroup rolls back the optimistic membership on error", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID), [makeGroup([])]);
    mockOnce(400, { error: "bad" });
    const { result } = renderHook(() => useAddIndexToGroup(EXPOSURE_ID, GROUP_ID), { wrapper });
    act(() => { result.current.mutate(101); });
    await waitFor(() => {
      const groups = client.getQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID)) ?? [];
      expect(groups[0].members).toEqual([]);
    });
  });

  // ---------------------- useRemoveIndexFromGroup ----------------------

  it("useRemoveIndexFromGroup filters indexId out of the group's members optimistically", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID), [makeGroup([101, 102])]);
    mockNever();
    const { result } = renderHook(() => useRemoveIndexFromGroup(EXPOSURE_ID, GROUP_ID), { wrapper });
    act(() => { result.current.mutate(101); });
    await waitFor(() => {
      const groups = client.getQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID)) ?? [];
      expect(groups[0].members).toEqual([102]);
    });
  });

  it("useRemoveIndexFromGroup replaces the group with the server-returned GroupEntry", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID), [makeGroup([101, 102])]);
    mockOnce(200, { ...makeGroup([102]) });
    const { result } = renderHook(() => useRemoveIndexFromGroup(EXPOSURE_ID, GROUP_ID), { wrapper });
    act(() => { result.current.mutate(101); });
    await waitFor(() => {
      const groups = client.getQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID)) ?? [];
      expect(groups[0].members).toEqual([102]);
    });
  });

  it("useRemoveIndexFromGroup rolls back the optimistic removal on error", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID), [makeGroup([101, 102])]);
    mockOnce(400, { error: "bad" });
    const { result } = renderHook(() => useRemoveIndexFromGroup(EXPOSURE_ID, GROUP_ID), { wrapper });
    act(() => { result.current.mutate(101); });
    await waitFor(() => {
      const groups = client.getQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID)) ?? [];
      expect(groups[0].members).toEqual([101, 102]);
    });
  });

  // ---------------------- useDeleteIndex ----------------------

  it("useDeleteIndex filters the index from BOTH indices and group members optimistically", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID),
      [makeIndex(101), makeIndex(102)]);
    client.setQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID),
      [makeGroup([101, 102])]);
    mockNever();
    const { result } = renderHook(() => useDeleteIndex(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(101); });
    await waitFor(() => {
      const indices = client.getQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID)) ?? [];
      const groups = client.getQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID)) ?? [];
      expect(indices.map((i) => i.id)).toEqual([102]);
      expect(groups[0].members).toEqual([102]);
    });
  });

  it("useDeleteIndex success leaves the optimistic state in place", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID),
      [makeIndex(101), makeIndex(102)]);
    client.setQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID),
      [makeGroup([101, 102])]);
    mockOnce(200, { deleted: 101 });
    const { result } = renderHook(() => useDeleteIndex(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(101); });
    await waitFor(() => expect(result.current.isPending).toBe(false));
    const indices = client.getQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID)) ?? [];
    const groups = client.getQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID)) ?? [];
    expect(indices.map((i) => i.id)).toEqual([102]);
    expect(groups[0].members).toEqual([102]);
  });

  it("useDeleteIndex rolls back BOTH indices and group caches on error", async () => {
    const { client, wrapper } = withClient();
    client.setQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID),
      [makeIndex(101), makeIndex(102)]);
    client.setQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID),
      [makeGroup([101, 102])]);
    mockOnce(403, { error: "only speculative indices can be deleted" });
    const { result } = renderHook(() => useDeleteIndex(EXPOSURE_ID), { wrapper });
    act(() => { result.current.mutate(101); });
    await waitFor(() => {
      const indices = client.getQueryData<IndexEntry[]>(queryKeys.indices(EXPOSURE_ID)) ?? [];
      const groups = client.getQueryData<GroupEntry[]>(queryKeys.groups(EXPOSURE_ID)) ?? [];
      expect(indices.map((i) => i.id)).toEqual([101, 102]);
      expect(groups[0].members).toEqual([101, 102]);
    });
  });
});
