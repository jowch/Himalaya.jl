import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type {
  SampleTagPair,
  PickerSampleRow,
  Sample,
  SampleTag,
  Trace,
  IndexEntry,
} from "../../src/api";

// ── navigate spy ─────────────────────────────────────────────────────────────
const navigateSpy = vi.fn();
vi.mock("react-router-dom", async (importOriginal) => {
  const actual = await importOriginal<typeof import("react-router-dom")>();
  return { ...actual, useNavigate: () => navigateSpy };
});

// ── helpers to build strict-typed fixtures ───────────────────────────────────
function tag(key: string, value: string): SampleTag {
  return { id: 1, key, value, source: "manual" };
}
function sample(id: number, name: string, tags: SampleTag[]): Sample {
  return { id, experiment_id: 1, name, display_name: name, notes: "", tags };
}
function pickerRow(s: Sample, exposureId: number): PickerSampleRow {
  return { sample: s, indexing_exposure_id: exposureId, all_exposures: [] };
}
function trace(): Trace {
  return { q: [0.1, 0.2], I: [1, 2], sigma: [0, 0] };
}
function indexEntry(exposureId: number, phase: string, score: number): IndexEntry {
  return {
    id: exposureId,
    exposure_id: exposureId,
    phase,
    basis: 1,
    score,
    r_squared: 0.9,
    lattice_d: null,
    ngc: null,
    status: "candidate",
    kind: "auto",
    inputs_hash: null,
    peaks: [],
    predicted_q: [],
  };
}

// ── mock data plane (mutated per test in beforeEach / inline) ─────────────────
const mutate = vi.fn();
let scopeState: { mutate: typeof mutate; isSuccess: boolean; error: Error | null } = {
  mutate,
  isSuccess: false,
  error: null,
};
let tagsState: { data: SampleTagPair[]; isLoading: boolean; isError: boolean };
let pickerState: { data: PickerSampleRow[]; isLoading: boolean; isError: boolean };

vi.mock("../../src/queries", () => ({
  useCorpusSampleTags: () => tagsState,
  useCorpusPickerSamples: () => pickerState,
  useScopeSeries: () => scopeState,
  useMemberTraces: (ids: number[]) =>
    new Map<number, Trace>(ids.map((id) => [id, trace()])),
  useMemberIndices: (ids: number[]) =>
    new Map<number, IndexEntry[]>(ids.map((id) => [id, [indexEntry(id, "Pn3m", 0.9)]])),
}));

// boneyard Skeleton: render children when not loading (mirror the folio harness).
vi.mock("boneyard-js/react", () => ({
  Skeleton: ({
    children,
    loading,
    fallback,
  }: {
    children: React.ReactNode;
    loading: boolean;
    fallback?: React.ReactNode;
  }) => (loading ? <>{fallback}</> : <>{children}</>),
}));

import { SeriesScopingPage } from "../../src/print/pages/SeriesScopingPage";

function renderPage(): void {
  const qc = new QueryClient();
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/series/new"]}>
        <SeriesScopingPage />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

function seed(): void {
  tagsState = {
    data: [{ key: "ratio", value: "1 : 0.5" }],
    isLoading: false,
    isError: false,
  };
  pickerState = {
    data: [
      pickerRow(sample(1, "A", [tag("ratio", "1 : 0")]), 37),
      pickerRow(sample(2, "B", [tag("ratio", "1 : 0.5")]), 65),
    ],
    isLoading: false,
    isError: false,
  };
}

beforeEach(() => {
  vi.clearAllMocks();
  scopeState = { mutate, isSuccess: false, error: null };
  seed();
});

describe("SeriesScopingPage", () => {
  it("renders a scope-sample-row per member", () => {
    renderPage();
    expect(screen.getAllByTestId("scope-sample-row")).toHaveLength(2);
  });

  it("surfaces a corpus-load error distinctly (not as an empty state)", () => {
    tagsState = { data: [], isLoading: false, isError: true };
    renderPage();
    expect(screen.getByText(/couldn't load the corpus/i)).toBeInTheDocument();
    expect(screen.queryByTestId("scope-sample-row")).not.toBeInTheDocument();
  });

  it("writes only confirmed members and navigates on build success", () => {
    renderPage();
    fireEvent.click(screen.getByRole("button", { name: /confirm & build/i }));
    expect(mutate).toHaveBeenCalledWith(
      expect.objectContaining({
        key: "ratio",
        tags: expect.arrayContaining([{ sampleId: 2, value: "1 : 0.5" }]),
      }),
    );
  });

  it("shows the ready foot line when no member is flagged", () => {
    renderPage();
    expect(screen.getByText(/All 2 values confirmed — ready to build/i)).toBeInTheDocument();
  });

  it("treats a sample missing the ordering key as a loose candidate, not a member", () => {
    // Sample A has no "ratio" value → it splits out as a loose match (a
    // candidate), leaving one clean member B. Build stays enabled.
    pickerState = {
      data: [
        pickerRow(sample(1, "A", []), 37),
        pickerRow(sample(2, "B", [tag("ratio", "1 : 0.5")]), 65),
      ],
      isLoading: false,
      isError: false,
    };
    renderPage();
    expect(screen.getAllByTestId("scope-sample-row")).toHaveLength(1);
    expect(screen.getByTestId("scope-candidate-row")).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /confirm & build/i })).not.toBeDisabled();
  });

  it("flagging a member's value warns, disables build; unflagging re-enables", () => {
    renderPage();
    const build = screen.getByRole("button", { name: /confirm & build/i });
    expect(build).not.toBeDisabled();
    // Click the first member's value control → re-open (flag) it.
    const flagButtons = screen.getAllByTestId("flag-button");
    fireEvent.click(flagButtons[0]!);
    expect(screen.getByText(/1 value to check before you can build/i)).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /confirm & build/i })).toBeDisabled();
    // Click it again → resolve; build re-enables.
    fireEvent.click(screen.getAllByTestId("flag-button")[0]!);
    expect(screen.getByRole("button", { name: /confirm & build/i })).not.toBeDisabled();
  });

  it("shows the no-ordering-variable empty state with a contact-sheet CTA on a cold corpus", () => {
    tagsState = { data: [], isLoading: false, isError: false };
    pickerState = {
      data: [pickerRow(sample(1, "A", []), 37)],
      isLoading: false,
      isError: false,
    };
    renderPage();
    expect(screen.getByTestId("empty-state")).toBeInTheDocument();
    expect(screen.queryByTestId("scope-sample-row")).not.toBeInTheDocument();
    const cta = screen.getByRole("button", { name: /contact sheet/i });
    fireEvent.click(cta);
    expect(navigateSpy).toHaveBeenCalledWith("/samples");
  });

  it("shows the write-error banner when the batch write fails", () => {
    scopeState = { mutate, isSuccess: false, error: new Error("boom") };
    renderPage();
    const banner = screen.getByRole("alert");
    expect(banner).toBeInTheDocument();
    expect(banner.textContent).toMatch(/could not write the scoping tags/i);
  });
});
