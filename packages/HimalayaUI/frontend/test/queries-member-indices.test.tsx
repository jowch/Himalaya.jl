import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { useMemberIndices } from "../src/queries";

function wrap(client = makeClient()) {
  return ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
}

describe("useMemberIndices", () => {
  beforeEach(() => vi.restoreAllMocks());
  it("fetches indices for each exposure id and maps them by id", async () => {
    vi.spyOn(global, "fetch").mockImplementation((input) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url.endsWith("/api/exposures/7/indices"))
        return Promise.resolve(new Response(JSON.stringify([
          { id: 1, exposure_id: 7, phase: "Pn3m", basis: 1, score: 0.9, r_squared: null,
            lattice_d: null, ngc: null, status: "candidate", kind: "auto",
            inputs_hash: null, peaks: [], predicted_q: [] },
        ]), { status: 200, headers: { "Content-Type": "application/json" } }));
      return Promise.resolve(new Response("[]", { status: 200, headers: { "Content-Type": "application/json" } }));
    });
    const { result } = renderHook(() => useMemberIndices([7]), { wrapper: wrap() });
    await waitFor(() => expect(result.current.get(7)?.[0]?.phase).toBe("Pn3m"));
  });
});
