import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor, act } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useGroups, useAddIndexToGroup, useRemoveIndexFromGroup, queryKeys,
} from "../src/queries";

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

describe("queries — groups", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("useGroups fetches when exposureId is provided", async () => {
    mockOnce(200, [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [] }]);
    const { wrapper } = withClient();
    const { result } = renderHook(() => useGroups(42), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(result.current.data).toHaveLength(1);
  });

  it("useGroups disabled when exposureId undefined", () => {
    const spy = vi.spyOn(global, "fetch");
    const { wrapper } = withClient();
    const { result } = renderHook(() => useGroups(undefined), { wrapper });
    expect(result.current.fetchStatus).toBe("idle");
    expect(spy).not.toHaveBeenCalled();
  });

  it("useAddIndexToGroup invalidates groups for exposure", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(200, { id: 2, exposure_id: 42, kind: "custom", active: true, members: [10, 11] });
    const { result } = renderHook(() => useAddIndexToGroup(42, 2), { wrapper });
    await act(async () => { await result.current.mutateAsync(11); });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.groups(42) });
  });

  it("useRemoveIndexFromGroup invalidates groups for exposure", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(200, { id: 2, exposure_id: 42, kind: "custom", active: true, members: [10] });
    const { result } = renderHook(() => useRemoveIndexFromGroup(42, 2), { wrapper });
    await act(async () => { await result.current.mutateAsync(11); });
    expect(invalidate).toHaveBeenCalledWith({ queryKey: queryKeys.groups(42) });
  });
});
