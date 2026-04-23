import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, fireEvent, waitFor } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { StaleIndicesBanner } from "../src/components/StaleIndicesBanner";

beforeEach(() => { vi.restoreAllMocks(); });

function mockResp(url: string, status: number, body: unknown): void {
  vi.spyOn(global, "fetch").mockImplementation(async (input) => {
    const u = typeof input === "string" ? input : (input as Request).url;
    if (u.endsWith(url)) {
      return new Response(JSON.stringify(body), {
        status, headers: { "Content-Type": "application/json" },
      });
    }
    return new Response("not found", { status: 404 });
  });
}

describe("<StaleIndicesBanner>", () => {
  it("renders nothing when exposureId is undefined", () => {
    const { container } = renderWithProviders(
      <StaleIndicesBanner exposureId={undefined} />,
    );
    expect(container.textContent).toBe("");
  });

  it("renders nothing when no indices are stale", async () => {
    mockResp("/api/exposures/42/indices", 200, [
      { id: 1, exposure_id: 42, phase: "Pn3m", basis: 0.1, score: 1,
        r_squared: 0.99, lattice_d: 10, status: "candidate" },
    ]);
    const { container } = renderWithProviders(
      <StaleIndicesBanner exposureId={42} />,
    );
    await waitFor(() => expect(container.textContent).toBe(""));
  });

  it("renders a re-analyze button when any index is stale", async () => {
    mockResp("/api/exposures/42/indices", 200, [
      { id: 1, exposure_id: 42, phase: "Pn3m", basis: 0.1, score: 1,
        r_squared: 0.99, lattice_d: 10, status: "stale" },
    ]);
    renderWithProviders(<StaleIndicesBanner exposureId={42} />);
    await waitFor(() =>
      expect(screen.getByRole("button", { name: /re-analyze/i })).toBeInTheDocument(),
    );
  });

  it("clicking re-analyze calls POST /api/exposures/:id/analyze", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(async (input) => {
      const u = typeof input === "string" ? input : (input as Request).url;
      if (u.endsWith("/api/exposures/42/indices")) {
        return new Response(JSON.stringify([
          { id: 1, exposure_id: 42, phase: "Pn3m", basis: 0.1, score: 1,
            r_squared: 0.99, lattice_d: 10, status: "stale" },
        ]), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (u.endsWith("/api/exposures/42/analyze")) {
        return new Response(JSON.stringify({ id: 42, analyzed: true }),
          { status: 200, headers: { "Content-Type": "application/json" } });
      }
      return new Response("not found", { status: 404 });
    });
    renderWithProviders(<StaleIndicesBanner exposureId={42} />);
    const btn = await screen.findByRole("button", { name: /re-analyze/i });
    fireEvent.click(btn);
    await waitFor(() => {
      const urls = fetchSpy.mock.calls.map((c) =>
        typeof c[0] === "string" ? c[0] : (c[0] as Request).url);
      expect(urls).toContain("/api/exposures/42/analyze");
    });
  });
});
