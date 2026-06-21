// test/SamplesPage.keyboard.test.tsx
//
// T2.4: page-level cursor + Corpus keyboard map.
//
// Pins:
//   ↑/↓ move the sample cursor (clamped, not circular).
//   Alt+↑/↓ = reorderUp/reorderDown — NOT bound on Corpus; cursor must NOT move.
//   Enter = openFocus → navigate `/sample/:id` (FLAT path, spec §6.4).
//   `r` = representative — NOT bound on Corpus; navigate must NOT be called.
//
// Helper names follow the project's existing RTL + vitest import conventions
// (see test/SamplesPage.staleSelect.test.tsx + test/ExperimentsHomePage.test.tsx).
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { MemoryRouter } from "react-router-dom";
import { SamplesPage } from "../src/print/pages/SamplesPage";
import * as queries from "../src/queries";

// ── navigate spy ──────────────────────────────────────────────────────────────
// Module-wide mock so every render in this file gets the same spy.
const navigate = vi.fn();
vi.mock("react-router-dom", async (orig) => {
  const m = await orig<typeof import("react-router-dom")>();
  return { ...m, useNavigate: () => navigate };
});

// ── Fixtures ──────────────────────────────────────────────────────────────────
// Three corpus samples (distinct ids so the cursor can step 0→1→2).
// No display_name — toSampleRowModel falls back to `name`.
const mkSample = (id: number) => ({
  id,
  experiment_id: 1,
  name: `Sample ${id}`,
  display_name: null as null,
  notes: null as null,
  tags: [] as [],
  q_units: "A-1",
});
const S0 = mkSample(10);
const S1 = mkSample(11);
const S2 = mkSample(12);

// ── Setup helper ──────────────────────────────────────────────────────────────
function renderSamplesPage() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });

  vi.spyOn(queries, "useCorpusSamples").mockReturnValue({
    data: [S0, S1, S2],
    isLoading: false,
    isError: false,
  } as unknown as ReturnType<typeof queries.useCorpusSamples>);

  // Empty exposures for all samples — keeps ThumbnailGallery from firing image
  // fetches (which fail under JSDOM's relative-URL resolver). The keyboard
  // cursor tests are sample-grain only; frame-level cursor tests don't require
  // real rendered thumbnails.
  const byId = new Map([
    [S0.id, [] as import("../src/api").Exposure[]],
    [S1.id, [] as import("../src/api").Exposure[]],
    [S2.id, [] as import("../src/api").Exposure[]],
  ]);
  vi.spyOn(queries, "useCorpusExposures").mockReturnValue({
    byId,
    isLoading: false,
  });
  vi.spyOn(queries, "useScreenedProgress").mockReturnValue({ screened: 0, total: 3 });
  vi.spyOn(queries, "useExperiments").mockReturnValue({
    data: [],
    isLoading: false,
  } as unknown as ReturnType<typeof queries.useExperiments>);
  vi.spyOn(queries, "useSetExposureStatusBatch").mockReturnValue({
    mutate: vi.fn(),
  } as unknown as ReturnType<typeof queries.useSetExposureStatusBatch>);

  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/samples"]}>
        <SamplesPage />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("SamplesPage — T2.4 page cursor + keyboard map", () => {
  beforeEach(() => {
    navigate.mockClear();
    vi.restoreAllMocks();
    navigate.mockClear(); // restore re-mocks module; re-clear after restore
  });

  it("ArrowDown moves the sample cursor forward; ArrowUp moves it back", () => {
    renderSamplesPage();

    // Cursor starts at sampleIndex=0 (S0, id=10).
    // ArrowDown → sampleIndex=1 (S1, id=11).
    fireEvent.keyDown(window, { key: "ArrowDown" });
    // ArrowDown again → sampleIndex=2 (S2, id=12).
    fireEvent.keyDown(window, { key: "ArrowDown" });
    // Enter → navigate to the ACTIVE sample (S2, id=12).
    fireEvent.keyDown(window, { key: "Enter" });
    expect(navigate).toHaveBeenCalledWith("/sample/12");

    navigate.mockClear();

    // ArrowUp → sampleIndex=1 (S1, id=11).
    fireEvent.keyDown(window, { key: "ArrowUp" });
    fireEvent.keyDown(window, { key: "Enter" });
    expect(navigate).toHaveBeenCalledWith("/sample/11");
  });

  it("Alt+ArrowUp/Down (reorderUp/reorderDown) do NOT move the sample cursor", () => {
    renderSamplesPage();

    // ArrowDown once → cursor at sampleIndex=1 (S1, id=11).
    fireEvent.keyDown(window, { key: "ArrowDown" });

    // Alt+ArrowUp = reorderUp — a DISTINCT combo, not bound on Corpus.
    fireEvent.keyDown(window, { key: "ArrowUp", altKey: true });
    // Alt+ArrowDown = reorderDown — same: not bound on Corpus.
    fireEvent.keyDown(window, { key: "ArrowDown", altKey: true });

    // Cursor should still be at sampleIndex=1 (S1, id=11).
    fireEvent.keyDown(window, { key: "Enter" });
    expect(navigate).toHaveBeenCalledWith("/sample/11");
  });

  it("ArrowDown clamps at the last sample", () => {
    renderSamplesPage();

    // Walk past the end.
    fireEvent.keyDown(window, { key: "ArrowDown" }); // → 1
    fireEvent.keyDown(window, { key: "ArrowDown" }); // → 2
    fireEvent.keyDown(window, { key: "ArrowDown" }); // stays 2 (clamped)

    fireEvent.keyDown(window, { key: "Enter" });
    expect(navigate).toHaveBeenCalledWith("/sample/12");
  });

  it("ArrowUp clamps at the first sample", () => {
    renderSamplesPage();

    // Already at 0; ArrowUp should not go negative.
    fireEvent.keyDown(window, { key: "ArrowUp" });

    fireEvent.keyDown(window, { key: "Enter" });
    expect(navigate).toHaveBeenCalledWith("/sample/10");
  });

  it("Enter (openFocus) navigates to /sample/:id (flat)", () => {
    renderSamplesPage();
    // Cursor at index 0 → sample id 10.
    fireEvent.keyDown(window, { key: "Enter" });
    expect(navigate).toHaveBeenCalledWith("/sample/10");
  });

  it("r key (representative) does nothing on Corpus — navigate not called", () => {
    renderSamplesPage();
    fireEvent.keyDown(window, { key: "r" });
    expect(navigate).not.toHaveBeenCalled();
  });

  it("l key (openLoupe) navigates to /sample/:id/loupe", () => {
    renderSamplesPage();
    // Cursor at index 0 → sample id 10.
    fireEvent.keyDown(window, { key: "l" });
    expect(navigate).toHaveBeenCalledWith("/sample/10/loupe");
  });
});
