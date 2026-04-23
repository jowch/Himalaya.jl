import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, waitFor, fireEvent } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { PhasePanel } from "../src/components/PhasePanel";
import { useAppState } from "../src/state";

beforeEach(() => { vi.restoreAllMocks(); });

function mockAll(indices: unknown[], groups: unknown[]): void {
  vi.spyOn(global, "fetch").mockImplementation(async (input) => {
    const u = typeof input === "string" ? input : (input as Request).url;
    if (u.endsWith("/indices")) {
      return new Response(JSON.stringify(indices),
        { status: 200, headers: { "Content-Type": "application/json" } });
    }
    if (u.endsWith("/groups")) {
      return new Response(JSON.stringify(groups),
        { status: 200, headers: { "Content-Type": "application/json" } });
    }
    if (u.match(/\/api\/groups\/\d+\/members\/\d+$/)) {
      return new Response(JSON.stringify({
        id: 2, exposure_id: 42, kind: "custom", active: true, members: [],
      }), { status: 200, headers: { "Content-Type": "application/json" } });
    }
    return new Response("not found", { status: 404 });
  });
}

describe("<PhasePanel> — active group", () => {
  it("renders a hint when no exposure is active", () => {
    renderWithProviders(<PhasePanel exposureId={undefined} />);
    expect(screen.getByText(/no exposure selected/i)).toBeInTheDocument();
  });

  it("renders each active-group index with phase, lattice, R² and a remove button", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 1.0,
          r_squared: 0.998, lattice_d: 12.5, status: "candidate",
          predicted_q: [0.7, 0.9], peaks: [] },
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
          r_squared: 0.71, lattice_d: 9.1, status: "candidate",
          predicted_q: [0.4, 0.6], peaks: [] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    await waitFor(() => expect(screen.getByText("Pn3m")).toBeInTheDocument());
    expect(screen.getByText(/12\.50/)).toBeInTheDocument();
    expect(screen.getByText(/0\.998/)).toBeInTheDocument();
    const rm = screen.getByRole("button", { name: /remove index 10/i });
    expect(rm).toBeInTheDocument();
  });

  it("clicking remove calls DELETE /api/groups/:gid/members/:indexId", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 1.0,
          r_squared: 0.998, lattice_d: 12.5, status: "candidate",
          predicted_q: [0.7, 0.9], peaks: [] },
      ],
      [{ id: 2, exposure_id: 42, kind: "custom", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const rm = await screen.findByRole("button", { name: /remove index 10/i });
    fireEvent.click(rm);
    await waitFor(() => {
      const spy = global.fetch as unknown as { mock: { calls: unknown[][] } };
      const urls = spy.mock.calls.map((c) =>
        typeof c[0] === "string" ? c[0] : (c[0] as Request).url);
      expect(urls.some((u) => u === "/api/groups/2/members/10")).toBe(true);
    });
  });
});

describe("<PhasePanel> — alternatives", () => {
  it("renders alternative indices with a + button", async () => {
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 1.0,
          r_squared: 0.998, lattice_d: 12.5, status: "candidate",
          predicted_q: [0.7, 0.9], peaks: [] },
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
          r_squared: 0.71, lattice_d: 9.1, status: "candidate",
          predicted_q: [0.4, 0.6], peaks: [] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    await waitFor(() => expect(screen.getByText("Im3m")).toBeInTheDocument());
    expect(screen.getByRole("button", { name: /add index 11/i })).toBeInTheDocument();
  });

  it("hovering an alternative sets hoveredIndexId; leaving clears it", async () => {
    useAppState.setState({ hoveredIndexId: undefined });
    mockAll(
      [
        { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 1.0,
          r_squared: 0.998, lattice_d: 12.5, status: "candidate",
          predicted_q: [0.7], peaks: [] },
        { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
          r_squared: 0.71, lattice_d: 9.1, status: "candidate",
          predicted_q: [0.4], peaks: [] },
      ],
      [{ id: 1, exposure_id: 42, kind: "auto", active: true, members: [10] }],
    );
    renderWithProviders(<PhasePanel exposureId={42} />);
    const row = await waitFor(() =>
      document.querySelector<HTMLElement>('[data-alternative-id="11"]') ?? (() => { throw new Error("not found"); })(),
    );
    fireEvent.mouseEnter(row);
    expect(useAppState.getState().hoveredIndexId).toBe(11);
    fireEvent.mouseLeave(row);
    expect(useAppState.getState().hoveredIndexId).toBeUndefined();
  });

  it("clicking + on an alternative posts to /api/groups/:gid/members", async () => {
    vi.spyOn(global, "fetch").mockImplementation(async (input) => {
      const u = typeof input === "string" ? input : (input as Request).url;
      if (u.endsWith("/indices")) {
        return new Response(JSON.stringify([
          { id: 10, exposure_id: 42, phase: "Pn3m", basis: 0.5, score: 1.0,
            r_squared: 0.998, lattice_d: 12.5, status: "candidate",
            predicted_q: [0.7], peaks: [] },
          { id: 11, exposure_id: 42, phase: "Im3m", basis: 0.3, score: 0.6,
            r_squared: 0.71, lattice_d: 9.1, status: "candidate",
            predicted_q: [0.4], peaks: [] },
        ]), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (u.endsWith("/groups")) {
        return new Response(JSON.stringify([
          { id: 2, exposure_id: 42, kind: "custom", active: true, members: [10] },
        ]), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (u.endsWith("/api/groups/2/members")) {
        return new Response(JSON.stringify({
          id: 2, exposure_id: 42, kind: "custom", active: true, members: [10, 11],
        }), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      return new Response("not found", { status: 404 });
    });

    renderWithProviders(<PhasePanel exposureId={42} />);
    const addBtn = await screen.findByRole("button", { name: /add index 11/i });
    fireEvent.click(addBtn);
    await waitFor(() => {
      const spy = global.fetch as unknown as { mock: { calls: unknown[][] } };
      const urls = spy.mock.calls.map((c) =>
        typeof c[0] === "string" ? c[0] : (c[0] as Request).url);
      expect(urls).toContain("/api/groups/2/members");
    });
  });
});
