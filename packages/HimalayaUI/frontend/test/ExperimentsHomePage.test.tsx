// test/ExperimentsHomePage.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ExperimentsHomePage } from "../src/pages/ExperimentsHomePage";
import * as api from "../src/api";

const navigate = vi.fn();
vi.mock("react-router-dom", async (orig) => {
  const m = await orig<typeof import("react-router-dom")>();
  return { ...m, useNavigate: () => navigate };
});

function renderPage() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter><ExperimentsHomePage /></MemoryRouter>
    </QueryClientProvider>,
  );
}

const EXP = (over: Partial<api.Experiment>): api.Experiment => ({
  id: 7, name: "Beamtime 1 · setup1", description: null, path: "/d",
  created_at: "2026-04-12T00:00:00",
  last_scanned_at: null, ingest_status: "complete",
  stats: { loads: 13, samples: 170, exposures: 682, sessions: 3, span_hours: 6.5, started_at: "2026-04-01T10:00:00" },
  ...over,
});

describe("ExperimentsHomePage (Phase E1)", () => {
  beforeEach(() => { navigate.mockClear(); });
  afterEach(() => vi.restoreAllMocks());

  it("lists experiment cards", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([EXP({}), EXP({ id: 8, name: "JC plate" })]);
    renderPage();
    expect(await screen.findByText("Beamtime 1 · setup1")).toBeInTheDocument();
    expect(screen.getByText("JC plate")).toBeInTheDocument();
  });

  it("shows when each beamtime collected its samples and how many are left to review", async () => {
    const timeline: api.CollectionTimeline = { start: "2026-04-01T10:00:00", samples: [[0, 60, 1], [120, 180, 0]], beam_off: [] };
    vi.spyOn(api, "listExperiments").mockResolvedValue([EXP({ stats: { ...EXP({}).stats!, review: { reviewed: 1, total: 2 }, timeline } })]);
    renderPage();
    const card = await screen.findByTestId("experiment-card-7");
    expect(card.querySelector("[data-testid='collection-ribbon'] [role='img']")).toHaveAccessibleName(/1 of 2 samples reviewed/);
    expect(screen.getByText("1 to review")).toBeInTheDocument();
  });

  it("opens the experiment on card click", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([EXP({})]);
    renderPage();
    fireEvent.click(await screen.findByTestId("experiment-card-7"));
    expect(navigate).toHaveBeenCalledWith("/experiments/7");
  });

  it("routes to /experiments/new on the New CTA", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([EXP({})]);
    renderPage();
    fireEvent.click(await screen.findByTestId("new-experiment-cta"));
    expect(navigate).toHaveBeenCalledWith("/experiments/new");
  });

  it("keeps the header New CTA working on the empty home page", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
    renderPage();
    await screen.findByText(/no experiments yet/i);
    const cta = screen.getByTestId("new-experiment-cta");
    expect(cta).toBeEnabled();
    fireEvent.click(cta);
    expect(navigate).toHaveBeenCalledWith("/experiments/new");
  });

  it("shows an empty state when there are no experiments", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
    renderPage();
    expect(await screen.findByText(/no experiments yet/i)).toBeInTheDocument();
  });

  it("renders the counts from stats", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({ stats: { loads: 13, samples: 170, exposures: 682, sessions: 3, span_hours: 6.5, started_at: "2026-04-01T10:00:00" } }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-7");
    expect(screen.getByText("170")).toBeInTheDocument();
    expect(screen.getByText("682")).toBeInTheDocument();
    expect(screen.getByText("13")).toBeInTheDocument();
  });

  it("shows the Up to date chip when complete", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({ ingest_status: "complete" }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-7");
    expect(screen.getByText("Up to date")).toBeInTheDocument();
  });

  it("shows a warning chip when the last ingest failed", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({ ingest_status: "failed" }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-7");
    expect(screen.getByText("Last ingest failed")).toBeInTheDocument();
  });

  it("shows the Scanning… chip when ingest_status is scanning", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({ ingest_status: "scanning" }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-7");
    expect(screen.getByText("Scanning…")).toBeInTheDocument();
  });

  it("groups experiments under year headings", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({
        id: 1, name: "Exp 2026",
        stats: { loads: 1, samples: 1, exposures: 1, sessions: 1, span_hours: 0, started_at: "2026-01-01T00:00:00" },
      }),
      EXP({
        id: 2, name: "Exp 2025",
        stats: { loads: 1, samples: 1, exposures: 1, sessions: 1, span_hours: 0, started_at: "2025-06-01T00:00:00" },
      }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-1");
    // Year headings are h2 elements; both 2026 and 2025 appear in the rail AND
    // as headings. Confirm the heading role specifically.
    expect(screen.getByRole("heading", { level: 2, name: "2026" })).toBeInTheDocument();
    expect(screen.getByRole("heading", { level: 2, name: "2025" })).toBeInTheDocument();
  });

  it("falls back to created_at year when stats.started_at is null", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({
        id: 5, name: "FallbackYear",
        created_at: "2024-03-15T00:00:00",
        stats: { loads: 0, samples: 0, exposures: 0, sessions: 0, span_hours: 0, started_at: null },
      }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-5");
    expect(screen.getByRole("heading", { level: 2, name: "2024" })).toBeInTheDocument();
  });

  it("shows 'indexing…' instead of counts when scanning", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({ ingest_status: "scanning" }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-7");
    expect(screen.getByText("indexing…")).toBeInTheDocument();
  });

  it("moves between experiment cards with the arrow keys", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([EXP({}), EXP({ id: 8, name: "JC plate" }), EXP({ id: 9, name: "Third" })]);
    renderPage();
    const first = await screen.findByTestId("experiment-card-7");
    const second = screen.getByTestId("experiment-card-8");
    const third = screen.getByTestId("experiment-card-9");
    fireEvent.keyDown(document.body, { key: "ArrowDown" });
    expect(document.activeElement).toBe(first);
    fireEvent.keyDown(first, { key: "ArrowRight" });
    expect(document.activeElement).toBe(second);
    fireEvent.keyDown(second, { key: "ArrowDown" });
    expect(document.activeElement).toBe(third);
    fireEvent.keyDown(third, { key: "ArrowDown" });
    expect(document.activeElement).toBe(third);
    fireEvent.keyDown(third, { key: "ArrowUp" });
    expect(document.activeElement).toBe(second);
    fireEvent.keyDown(second, { key: "ArrowLeft" });
    expect(document.activeElement).toBe(first);
  });

  it("leaves arrow keys alone while a dialog is open", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([EXP({})]);
    renderPage();
    await screen.findByTestId("experiment-card-7");
    const dialog = document.createElement("div");
    dialog.setAttribute("role", "dialog");
    dialog.setAttribute("aria-modal", "true");
    document.body.appendChild(dialog);
    fireEvent.keyDown(document.body, { key: "ArrowDown" });
    expect(document.activeElement).toBe(document.body);
    dialog.remove();
  });
});
