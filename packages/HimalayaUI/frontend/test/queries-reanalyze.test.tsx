import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { useReanalyzeExposure } from "../src/queries";
import { useQueueOpStatus } from "../src/lib/queue/hooks";
import { pendingDeferreds } from "../src/lib/queue/deferred";

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

function mockNever(): void {
  vi.spyOn(global, "fetch").mockImplementation(() => new Promise(() => {}));
}

describe("queries — useReanalyzeExposure (queue-driven, M2.5)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
  });

  it("calls api.reanalyzeExposure with the exposure id", async () => {
    const { wrapper } = withClient();
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValueOnce(
      new Response(JSON.stringify({ id: EXPOSURE_ID, analyzed: true }), {
        status: 200, headers: { "Content-Type": "application/json" },
      }),
    );
    const { result } = renderHook(
      () => useReanalyzeExposure(EXPOSURE_ID), { wrapper },
    );
    act(() => { result.current.mutate({}); });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(fetchSpy).toHaveBeenCalled();
    const url = (fetchSpy.mock.calls[0]?.[0] ?? "").toString();
    expect(url).toBe(`/api/exposures/${EXPOSURE_ID}/analyze`);
  });

  it("isPending is true while the request is in flight", async () => {
    const { wrapper } = withClient();
    mockNever();
    const { result } = renderHook(
      () => useReanalyzeExposure(EXPOSURE_ID), { wrapper },
    );
    act(() => { result.current.mutate({}); });
    await waitFor(() => expect(result.current.isPending).toBe(true));
  });

  it("isPending flips false after the request settles", async () => {
    const { wrapper } = withClient();
    mockOnce(200, { id: EXPOSURE_ID, analyzed: true });
    const { result } = renderHook(
      () => useReanalyzeExposure(EXPOSURE_ID), { wrapper },
    );
    act(() => { result.current.mutate({}); });
    await waitFor(() => expect(result.current.isPending).toBe(false));
    expect(result.current.isSuccess).toBe(true);
  });

  it("does not invalidate peaks/indices/groups on success (SSE post_state owns the cache update)", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(200, { id: EXPOSURE_ID, analyzed: true });
    const { result } = renderHook(
      () => useReanalyzeExposure(EXPOSURE_ID), { wrapper },
    );
    act(() => { result.current.mutate({}); });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(invalidate).not.toHaveBeenCalled();
  });

  it("useQueueOpStatus reflects the pending in-flight reanalyze", async () => {
    const { wrapper } = withClient();
    mockNever();
    const { result } = renderHook(
      () => {
        const m = useReanalyzeExposure(EXPOSURE_ID);
        const status = useQueueOpStatus("reanalyze_exposure", EXPOSURE_ID);
        return { m, status };
      },
      { wrapper },
    );
    expect(result.current.status).toBe("idle");
    act(() => { result.current.m.mutate({}); });
    await waitFor(() => expect(result.current.status).toBe("pending"));
  });
});
