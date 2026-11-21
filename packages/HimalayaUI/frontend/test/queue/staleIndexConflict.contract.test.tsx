import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { renderHook, act, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "../test-utils";
import { useAddAssignmentPhase, useRemoveAssignmentPhase, queryKeys } from "../../src/queries";
import { setToastImpl } from "../../src/lib/toast";
import { pendingDeferreds } from "../../src/lib/queue/deferred";
import type { Assignment } from "../../src/api";

const EXPOSURE_ID = 5;
const STALE_INDEX = 999;
const KEPT_INDEX = 7;
const BEFORE: Assignment = { exposure_id: EXPOSURE_ID, state: "indexed", members: [KEPT_INDEX] };
const CONFLICT_BODY = { error: "this index changed after re-analysis, so pick it again" };

describe("a phase call on a stale index (409 from POST /assignment/members)", () => {
  let toasts: string[] = [];
  let sent: Array<{ method: string; url: string; at: number }> = [];

  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
    toasts = [];
    sent = [];
    setToastImpl((msg) => { toasts.push(msg); });
  });
  afterEach(() => { setToastImpl(null); });

  function stubServer(addStatus: number): void {
    vi.spyOn(global, "fetch").mockImplementation(async (input, init) => {
      const url = String(input);
      const method = init?.method ?? "GET";
      sent.push({ method, url, at: Date.now() });
      if (method === "POST") {
        return new Response(JSON.stringify(addStatus === 409 ? CONFLICT_BODY : { error: "boom" }), {
          status: addStatus, headers: { "Content-Type": "application/json" },
        });
      }
      return new Response(JSON.stringify({ exposure_id: EXPOSURE_ID, state: "null", members: [] }), {
        status: 200, headers: { "Content-Type": "application/json" },
      });
    });
  }

  function render() {
    const client = makeClient();
    client.setQueryData(queryKeys.assignment(EXPOSURE_ID), BEFORE);
    const wrapper = ({ children }: { children: ReactNode }) => (
      <QueryClientProvider client={client}>{children}</QueryClientProvider>
    );
    const { result } = renderHook(
      () => ({ add: useAddAssignmentPhase(EXPOSURE_ID), remove: useRemoveAssignmentPhase(EXPOSURE_ID) }),
      { wrapper },
    );
    return { client, result };
  }

  it("rolls the call back once, says why, and lets the next edit through at once", async () => {
    stubServer(409);
    const { client, result } = render();
    act(() => {
      result.current.add.mutate(STALE_INDEX);
      result.current.remove.mutate(KEPT_INDEX);
    });
    await waitFor(() => expect(sent.some((r) => r.method === "DELETE")).toBe(true), { timeout: 2000 });
    const posts = sent.filter((r) => r.method === "POST");
    expect(posts).toHaveLength(1);
    expect(toasts).toContain(`Couldn't call phase: ${CONFLICT_BODY.error}`);
    await waitFor(() =>
      expect(client.getQueryData<Assignment>(queryKeys.assignment(EXPOSURE_ID))?.members).not.toContain(STALE_INDEX),
    );
  });

  it("a 500 on the same call is retried, so the edit behind it waits", async () => {
    stubServer(500);
    const { result } = render();
    act(() => {
      result.current.add.mutate(STALE_INDEX);
      result.current.remove.mutate(KEPT_INDEX);
    });
    await new Promise((r) => setTimeout(r, 300));
    expect(sent.some((r) => r.method === "DELETE")).toBe(false);
  });
});
