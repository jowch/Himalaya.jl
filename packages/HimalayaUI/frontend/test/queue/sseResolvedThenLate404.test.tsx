/**
 * SSE-resolves-first then late HTTP 404 — entire op must be suppressed.
 *
 * Race: a mutation is in flight. The server commits + broadcasts SSE before
 * the HTTP response is delivered (common under WAN latency or LAN buffering).
 * `replayCoordinator.handleRemoteEvent` resolves the deferred and aborts the
 * controller; later, a 404 from the server (or a 404 the server emitted for
 * a retry of a delete that already landed) hits the client.
 *
 * Expected end state:
 *   - Optimistic effect held (no rollback)
 *   - Mutation settled as SUCCESS via the SSE-synthesized response
 *   - No toast surfaced
 *   - pendingDeferreds registry empty
 *
 * Pre-existing `treats404AsSuccess.test.tsx` covers the flag in isolation
 * (404 alone keeps the optimistic state). This file covers the orthogonal
 * race ordering — SSE first, then 404 — that the flag is downstream of:
 * once `deferred.resolve(synth)` lands, NO 404 path runs, so the flag is
 * irrelevant to correctness here. Sibling to `rollbackSymmetry.test.ts` per
 * issue #125 item 3.
 */
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "../test-utils";
import { useRemovePeak, queryKeys } from "../../src/queries";
import { setToastImpl } from "../../src/lib/toast";
import { pendingDeferreds } from "../../src/lib/queue/deferred";
import { handleRemoteEvent } from "../../src/lib/queue/replayCoordinator";
import type { Peak } from "../../src/api";

const EXPOSURE_ID = 5;
const PEAK: Peak = {
  id: 7, exposure_id: EXPOSURE_ID, q: 0.5, intensity: 1.0,
  prominence: 0.8, sharpness: 30, source: "manual", excluded: false,
};

describe("SSE-resolved-then-late-404 race — useRemovePeak", () => {
  let toastCalls: Array<{ msg: string; kind: string }>;

  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
    toastCalls = [];
    setToastImpl((msg, kind) => { toastCalls.push({ msg, kind }); });
  });

  afterEach(() => {
    setToastImpl(null);
  });

  it("optimistic state held, mutation succeeds, no rollback, no toast", async () => {
    // Mock fetch under our control so we can deliver the 404 strictly AFTER
    // the SSE has won the race.
    let resolveFetch!: (response: Response) => void;
    const fetchPromise = new Promise<Response>((res) => { resolveFetch = res; });
    vi.spyOn(global, "fetch").mockImplementation(() => fetchPromise);

    const client = makeClient();
    client.setQueryData(queryKeys.peaks(EXPOSURE_ID), [PEAK]);
    const wrapper = ({ children }: { children: ReactNode }) => (
      <QueryClientProvider client={client}>{children}</QueryClientProvider>
    );
    const { result } = renderHook(() => useRemovePeak(EXPOSURE_ID), { wrapper });

    act(() => { result.current.mutate(PEAK.id); });

    // Optimistic effect — peak gone before either SSE or HTTP returns.
    expect(client.getQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID))).toEqual([]);

    // Wait for mutationFn to register its deferred. The clientOpId is minted
    // inside mutate() and isn't visible to the test, so we grab it off the
    // registry as a signal that the queue has wired up the race.
    await waitFor(() => expect(pendingDeferreds.size).toBe(1));
    const clientOpId = [...pendingDeferreds.keys()][0]!;

    // SSE wins. replayCoordinator resolves the deferred with a synthesized
    // response AND aborts the controller; clearDeferred drops the entry.
    handleRemoteEvent(
      {
        id: 11,
        kind: "peak_removed",
        entity_type: "exposure",
        entity_id: EXPOSURE_ID,
        client_op_id: clientOpId,
        payload: { peakId: PEAK.id, q: PEAK.q },
      } as any,
      client,
      client.getMutationCache(),
    );
    expect(pendingDeferreds.size).toBe(0);

    // The late 404 lands AFTER SSE. The mock fetch ignores AbortSignal — that
    // is the race we're pinning: the response was already in flight when the
    // controller aborted. The HTTP wire-up in useQueueMutation routes this
    // through `.catch(deferred.reject)`, but `deferred` already resolved, so
    // the reject is a no-op. onSuccess fires with the SSE-synth response;
    // onError never runs.
    resolveFetch(
      new Response(JSON.stringify({ error: "peak not found" }), {
        status: 404, headers: { "Content-Type": "application/json" },
      }),
    );

    await waitFor(() => expect(result.current.isPending).toBe(false));

    // Optimistic delete sticks. (treats404AsSuccess would also keep this,
    // so the remaining assertions distinguish the SSE-wins path.)
    expect(client.getQueryData<Peak[]>(queryKeys.peaks(EXPOSURE_ID))).toEqual([]);

    // The mutation settled as a SUCCESS — not as a swallowed 404 error.
    // Under treats404AsSuccess the hook would expose the ApiError on
    // `error` (or undefined depending on retry); under SSE-wins, success.
    expect(result.current.isSuccess).toBe(true);
    expect(result.current.error).toBeNull();

    // No toast: neither validation nor infrastructure surfaced.
    expect(toastCalls).toEqual([]);
  });
});
