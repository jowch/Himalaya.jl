import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { useAssignment, queryKeys } from "../src/queries";
import * as api from "../src/api";

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

describe("queries — assignment (Plan D-1)", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("queryKeys.assignment is exposure-scoped and distinct from indices", () => {
    expect(queryKeys.assignment(42)).toEqual(["exposure", 42, "assignment"]);
    expect(queryKeys.assignment(42)).not.toEqual(queryKeys.indices(42));
  });

  it("useAssignment fetches when exposureId is provided", async () => {
    mockOnce(200, { exposure_id: 42, state: "indexed", members: [10, 11] });
    const { wrapper } = withClient();
    const { result } = renderHook(() => useAssignment(42), { wrapper });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(result.current.data).toEqual({ exposure_id: 42, state: "indexed", members: [10, 11] });
  });

  it("useAssignment disabled when exposureId undefined", () => {
    const spy = vi.spyOn(global, "fetch");
    const { wrapper } = withClient();
    const { result } = renderHook(() => useAssignment(undefined), { wrapper });
    expect(result.current.fetchStatus).toBe("idle");
    expect(spy).not.toHaveBeenCalled();
  });
});

describe("api — assignment fetchers (Plan D-1)", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("getAssignment GETs the assignment endpoint", async () => {
    const spy = vi.spyOn(global, "fetch").mockResolvedValueOnce(
      new Response(JSON.stringify({ exposure_id: 7, state: "indexed", members: [] }),
        { status: 200, headers: { "Content-Type": "application/json" } }));
    const got = await api.getAssignment(7);
    expect(got.exposure_id).toBe(7);
    const [path, init] = spy.mock.calls[0]!;
    expect(path).toBe("/api/exposures/7/assignment");
    expect((init as RequestInit).method).toBe("GET");
  });

  it("addAssignmentPhase POSTs index_id to /assignment/members", async () => {
    const spy = vi.spyOn(global, "fetch").mockResolvedValueOnce(
      new Response(JSON.stringify({ exposure_id: 7, state: "indexed", members: [10] }),
        { status: 200, headers: { "Content-Type": "application/json" } }));
    await api.addAssignmentPhase(7, 10, { username: "alice" });
    const [path, init] = spy.mock.calls[0]!;
    expect(path).toBe("/api/exposures/7/assignment/members");
    expect((init as RequestInit).method).toBe("POST");
    expect(JSON.parse((init as RequestInit).body as string)).toEqual({ index_id: 10 });
    expect((init as RequestInit).headers).toMatchObject({ "X-Username": "alice" });
  });

  it("removeAssignmentPhase DELETEs /assignment/members/:index_id", async () => {
    const spy = vi.spyOn(global, "fetch").mockResolvedValueOnce(
      new Response(JSON.stringify({ exposure_id: 7, state: "indexed", members: [] }),
        { status: 200, headers: { "Content-Type": "application/json" } }));
    await api.removeAssignmentPhase(7, 10, { username: "alice" });
    const [path, init] = spy.mock.calls[0]!;
    expect(path).toBe("/api/exposures/7/assignment/members/10");
    expect((init as RequestInit).method).toBe("DELETE");
  });

  it("createCustomIndex POSTs {phase, basis, ratios} to /custom-index", async () => {
    const spy = vi.spyOn(global, "fetch").mockResolvedValueOnce(
      new Response(JSON.stringify({ id: 77, exposure_id: 7, phase: "Pn3m", basis: 0.15,
        score: null, r_squared: null, lattice_d: 197, ngc: null, status: "candidate",
        kind: "speculative", inputs_hash: "h", peaks: [], predicted_q: [0.15],
        event_id: 1, view_row_id: 1 }),
        { status: 200, headers: { "Content-Type": "application/json" } }));
    await api.createCustomIndex(7, "Pn3m", 0.15, [1, 1.2247448714], { username: "alice" });
    const [path, init] = spy.mock.calls[0]!;
    expect(path).toBe("/api/exposures/7/custom-index");
    expect((init as RequestInit).method).toBe("POST");
    // `ratios` is load-bearing, not decoration: without it the backend scans
    // the full core ratio series and claims reflections the modal never drew.
    // Ratios rather than a count because the modal's series is not a positional
    // prefix of the backend's for Hexagonal.
    expect(JSON.parse((init as RequestInit).body as string))
      .toEqual({ phase: "Pn3m", basis: 0.15, ratios: [1, 1.2247448714] });
  });

  it("setAssignmentState POSTs state to /assignment/state", async () => {
    const spy = vi.spyOn(global, "fetch").mockResolvedValueOnce(
      new Response(JSON.stringify({ exposure_id: 7, state: "form_factor", members: [] }),
        { status: 200, headers: { "Content-Type": "application/json" } }));
    await api.setAssignmentState(7, "form_factor", { username: "alice" });
    const [path, init] = spy.mock.calls[0]!;
    expect(path).toBe("/api/exposures/7/assignment/state");
    expect((init as RequestInit).method).toBe("POST");
    expect(JSON.parse((init as RequestInit).body as string)).toEqual({ state: "form_factor" });
  });
});
