import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, act } from "@testing-library/react";
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

function renderPage(seedSampleIds?: number[]): { rerender: () => void } {
  const qc = new QueryClient();
  const entry =
    seedSampleIds === undefined
      ? "/series/new"
      : { pathname: "/series/new", state: { seedSampleIds } };
  const tree = (): JSX.Element => (
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[entry]}>
        <SeriesScopingPage />
      </MemoryRouter>
    </QueryClientProvider>
  );
  const result = render(tree());
  // A FRESH element reference each rerender so React re-renders (and re-runs the
  // success effect after `scopeState.isSuccess` is flipped) rather than bailing
  // on an identical element reference.
  return { rerender: () => result.rerender(tree()) };
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

/** A 3-member corpus (all with a value) for skip/payload-exclusion cases. */
function seed3(): void {
  tagsState = { data: [{ key: "ratio", value: "1 : 0.5" }], isLoading: false, isError: false };
  pickerState = {
    data: [
      pickerRow(sample(1, "A", [tag("ratio", "1 : 0")]), 37),
      pickerRow(sample(2, "B", [tag("ratio", "1 : 0.5")]), 65),
      pickerRow(sample(3, "C", [tag("ratio", "1 : 1")]), 66),
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

  it("writes every member's read and navigates to /series on build success", () => {
    const { rerender } = renderPage();
    fireEvent.click(screen.getByRole("button", { name: /confirm & build/i }));
    expect(mutate).toHaveBeenCalledWith(
      expect.objectContaining({
        key: "ratio",
        tags: expect.arrayContaining([
          { sampleId: 1, value: "1 : 0" },
          { sampleId: 2, value: "1 : 0.5" },
        ]),
      }),
    );
    // No navigation until the write actually settles.
    expect(navigateSpy).not.toHaveBeenCalledWith("/series");
    // The write succeeds → the success effect navigates to the folio.
    scopeState = { mutate, isSuccess: true, error: null };
    act(() => rerender());
    expect(navigateSpy).toHaveBeenCalledWith("/series");
  });

  it("shows the ready foot line with the kept count when nothing is skipped", () => {
    renderPage();
    expect(screen.getByText(/2 values ready to commit/i)).toBeInTheDocument();
  });

  it("treats a sample missing the ordering key as a candidate, not a member", () => {
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
    expect(screen.getByTestId("scope-candidate")).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /confirm & build/i })).not.toBeDisabled();
  });

  it("skipping a member excludes it from the write but does not block the build", () => {
    seed3();
    renderPage();
    // Members sort low→high: A(1:0), B(1:0.5), C(1:1). Skip B (index 1).
    const flagButtons = screen.getAllByTestId("flag-button");
    fireEvent.click(flagButtons[1]!);
    // Build stays enabled (A and C are still kept).
    const build = screen.getByRole("button", { name: /confirm & build/i });
    expect(build).not.toBeDisabled();
    // Foot line annotates the skip.
    expect(screen.getByText(/2 values ready to commit · 1 skipped/i)).toBeInTheDocument();
    // The write excludes the skipped member B.
    fireEvent.click(build);
    expect(mutate).toHaveBeenCalledWith(
      expect.objectContaining({
        key: "ratio",
        tags: [
          { sampleId: 1, value: "1 : 0" },
          { sampleId: 3, value: "1 : 1" },
        ],
      }),
    );
  });

  it("disables build when every member is skipped, and re-enables on unskip", () => {
    renderPage();
    const flagButtons = screen.getAllByTestId("flag-button");
    fireEvent.click(flagButtons[0]!);
    fireEvent.click(screen.getAllByTestId("flag-button")[1]!);
    // Both members skipped → nothing to commit.
    expect(screen.getByText(/keep at least one value to build/i)).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /confirm & build/i })).toBeDisabled();
    // Unskip one → build re-enables.
    fireEvent.click(screen.getAllByTestId("flag-button")[0]!);
    expect(screen.getByRole("button", { name: /confirm & build/i })).not.toBeDisabled();
  });

  it("renders candidates as informational discovery with no add control", () => {
    pickerState = {
      data: [
        pickerRow(sample(1, "A", [tag("ratio", "1 : 0")]), 37),
        pickerRow(sample(2, "B", []), 65),
      ],
      isLoading: false,
      isError: false,
    };
    renderPage();
    // One member (A), one candidate (B).
    expect(screen.getAllByTestId("scope-sample-row")).toHaveLength(1);
    expect(screen.getByTestId("scope-candidate")).toBeInTheDocument();
    // No "+ Add to series" control anywhere.
    expect(screen.queryByRole("button", { name: /add to series/i })).toBeNull();
    expect(screen.queryByText(/add to series/i)).toBeNull();
    // The member count is unchanged (no way to fold a candidate in).
    expect(screen.getAllByTestId("scope-sample-row")).toHaveLength(1);
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

  it("scopes the proposal to the seeded sample ids when arrived with a seed", () => {
    // 3-member corpus (A,B,C); seed only A(1) and C(3). The proposal — and thus
    // the rendered members + the write — must scope to just the seeded samples.
    seed3();
    renderPage([1, 3]);
    expect(screen.getAllByTestId("scope-sample-row")).toHaveLength(2);
    fireEvent.click(screen.getByRole("button", { name: /confirm & build/i }));
    expect(mutate).toHaveBeenCalledWith(
      expect.objectContaining({
        key: "ratio",
        tags: [
          { sampleId: 1, value: "1 : 0" },
          { sampleId: 3, value: "1 : 1" },
        ],
      }),
    );
  });

  it("scopes to the whole corpus on a direct (unseeded) visit", () => {
    // Same 3-member corpus, no seed → all three are members.
    seed3();
    renderPage();
    expect(screen.getAllByTestId("scope-sample-row")).toHaveLength(3);
  });

  it("shows the write-error banner when the batch write fails", () => {
    scopeState = { mutate, isSuccess: false, error: new Error("boom") };
    renderPage();
    const banner = screen.getByRole("alert");
    expect(banner).toBeInTheDocument();
    expect(banner.textContent).toMatch(/could not write the scoping tags/i);
  });
});
