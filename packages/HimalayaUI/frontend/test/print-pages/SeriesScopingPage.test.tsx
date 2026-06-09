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
  Series,
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

// A minimal full `Series` (api.isFullSeries true) — the created series Op B
// (createSeries) resolves with; the page reads its id for navigation.
function fullSeries(id: number): Series {
  return {
    id,
    title: "Series by ratio",
    description: null,
    content_hash: "h",
    created_by: null,
    created_at: null,
    updated_at: null,
    forked_from_id: null,
    forked_at_hash: null,
    forked_from_title: null,
    view_grouping_mode: null,
    view_show_peak_ticks: null,
    view_show_peak_labels: null,
    ordering_variable: "ratio",
    order_rule: "ascending",
    state: "draft",
    members: [],
    samples: [],
  };
}

// ── mock data plane (mutated per test in beforeEach / inline) ─────────────────
// Two separate single-write queue ops (the M-A scope→create chain):
//   Op A  scopeSeries.mutate({key, tags})    — the tag batch write
//   Op B  createSeries.mutate({title, …})    — POST /api/series
type OpState = {
  mutate: ReturnType<typeof vi.fn>;
  isSuccess: boolean;
  error: Error | null;
  data: Series | undefined;
};
const scopeMutate = vi.fn();
const createMutate = vi.fn();
let scopeState: OpState = { mutate: scopeMutate, isSuccess: false, error: null, data: undefined };
let createState: OpState = { mutate: createMutate, isSuccess: false, error: null, data: undefined };
let tagsState: { data: SampleTagPair[]; isLoading: boolean; isError: boolean };
let pickerState: { data: PickerSampleRow[]; isLoading: boolean; isError: boolean };

vi.mock("../../src/queries", () => ({
  useCorpusSampleTags: () => tagsState,
  useCorpusPickerSamples: () => pickerState,
  useScopeSeries: () => scopeState,
  useCreateSeries: () => createState,
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
  scopeState = { mutate: scopeMutate, isSuccess: false, error: null, data: undefined };
  createState = { mutate: createMutate, isSuccess: false, error: null, data: undefined };
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

  it("warm: writes the tags (Op A), then creates the series (Op B), then navigates to /series/:id", () => {
    const { rerender } = renderPage();
    fireEvent.click(screen.getByRole("button", { name: /confirm & build/i }));
    // Op A fires FIRST — the tag batch write (no create body on this op).
    expect(scopeMutate).toHaveBeenCalledWith(
      expect.objectContaining({
        key: "ratio",
        tags: expect.arrayContaining([
          { sampleId: 1, value: "1 : 0" },
          { sampleId: 2, value: "1 : 0.5" },
        ]),
      }),
    );
    // Op B has NOT fired yet (gated on Op A's success).
    expect(createMutate).not.toHaveBeenCalled();
    expect(navigateSpy).not.toHaveBeenCalled();

    // Op A lands → the chain fires Op B (the create) with the scoped body.
    scopeState = { mutate: scopeMutate, isSuccess: true, error: null, data: undefined };
    act(() => rerender());
    expect(createMutate).toHaveBeenCalledWith(
      expect.objectContaining({
        title: "Series by ratio",
        ordering_variable: "ratio",
        samples: expect.arrayContaining([
          { sample_id: 1, position: 0 },
          { sample_id: 2, position: 1 },
        ]),
      }),
    );
    // Still no navigation until the create itself settles.
    expect(navigateSpy).not.toHaveBeenCalled();

    // Op B lands with the created Series → navigate to the new builder.
    createState = { mutate: createMutate, isSuccess: true, error: null, data: fullSeries(7) };
    act(() => rerender());
    expect(navigateSpy).toHaveBeenCalledWith("/series/7");
  });

  it("navigates to /series/:id even when Op B resolves with the SSE-win PARTIAL shape (id only, no members/state)", () => {
    // Under the live SSE-wins race, the queue resolves the create's deferred
    // with saveSeriesMutator.synthesizeFromSse → a PARTIAL series ({...base,
    // ...payload, id: entity_id}) — NO members/state, so isFullSeries is false,
    // but the new id IS present. The success effect must read the id directly,
    // not gate on isFullSeries (which would yield "/series", the folio). Mocked
    // tests normally drain SSE so HTTP wins and data is full — this asserts the
    // shape the full-series mock structurally hides.
    const { rerender } = renderPage();
    fireEvent.click(screen.getByRole("button", { name: /confirm & build/i }));
    // Op A lands → Op B fires.
    scopeState = { mutate: scopeMutate, isSuccess: true, error: null, data: undefined };
    act(() => rerender());
    expect(createMutate).toHaveBeenCalled();
    // Op B resolves with the PARTIAL SSE-synthesized shape: id only.
    createState = {
      mutate: createMutate,
      isSuccess: true,
      error: null,
      data: { id: 9 } as unknown as Series,
    };
    act(() => rerender());
    expect(navigateSpy).toHaveBeenCalledWith("/series/9");
    expect(navigateSpy).not.toHaveBeenCalledWith("/series");
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
    // The tag write (Op A) excludes the skipped member B.
    fireEvent.click(build);
    expect(scopeMutate).toHaveBeenCalledWith(
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

  it("cold path: names the variable, writes the tags (Op A), creates the series (Op B), navigates to /series/:id", () => {
    // No proposable key (empty tags) + a deliberate seed → the cold-assign
    // path. Naming the variable and assigning every value opens the build,
    // which runs the same two-op scope→create chain as the warm path.
    tagsState = { data: [], isLoading: false, isError: false };
    pickerState = {
      data: [
        pickerRow(sample(1, "A", []), 37),
        pickerRow(sample(2, "B", []), 65),
      ],
      isLoading: false,
      isError: false,
    };
    const { rerender } = renderPage([1, 2]);
    // Name the ordering variable.
    fireEvent.change(screen.getByTestId("cold-key-input"), {
      target: { value: "lipid ratio" },
    });
    // Assign every sample's value.
    const valueInputs = screen.getAllByTestId("cold-assign-row");
    fireEvent.change(valueInputs[0]!.querySelector("input")!, { target: { value: "1:0" } });
    fireEvent.change(valueInputs[1]!.querySelector("input")!, { target: { value: "1:1" } });
    const build = screen.getByRole("button", { name: /confirm & build/i });
    expect(build).not.toBeDisabled();
    fireEvent.click(build);
    // Op A — the tag batch write under the named key.
    expect(scopeMutate).toHaveBeenCalledWith(
      expect.objectContaining({
        key: "lipid ratio",
        tags: [
          { sampleId: 1, value: "1:0" },
          { sampleId: 2, value: "1:1" },
        ],
      }),
    );
    expect(createMutate).not.toHaveBeenCalled();

    // Op A lands → Op B fires with the create body built from the cold rows.
    scopeState = { mutate: scopeMutate, isSuccess: true, error: null, data: undefined };
    act(() => rerender());
    expect(createMutate).toHaveBeenCalledWith(
      expect.objectContaining({
        title: "Series by lipid ratio",
        ordering_variable: "lipid ratio",
        samples: [
          { sample_id: 1, position: 0 },
          { sample_id: 2, position: 1 },
        ],
      }),
    );

    // Op B lands with the created Series → navigate to the new builder.
    createState = { mutate: createMutate, isSuccess: true, error: null, data: fullSeries(9) };
    act(() => rerender());
    expect(navigateSpy).toHaveBeenCalledWith("/series/9");
  });

  it("scopes the proposal to the seeded sample ids when arrived with a seed", () => {
    // 3-member corpus (A,B,C); seed only A(1) and C(3). The proposal — and thus
    // the rendered members + the write — must scope to just the seeded samples.
    seed3();
    renderPage([1, 3]);
    expect(screen.getAllByTestId("scope-sample-row")).toHaveLength(2);
    fireEvent.click(screen.getByRole("button", { name: /confirm & build/i }));
    expect(scopeMutate).toHaveBeenCalledWith(
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

  it("shows the error banner when the tag write (Op A) fails", () => {
    scopeState = { mutate: scopeMutate, isSuccess: false, error: new Error("boom"), data: undefined };
    renderPage();
    const banner = screen.getByRole("alert");
    expect(banner).toBeInTheDocument();
    expect(banner.textContent).toMatch(/could not build the series/i);
  });

  it("shows the error banner when the series create (Op B) fails", () => {
    createState = { mutate: createMutate, isSuccess: false, error: new Error("boom"), data: undefined };
    renderPage();
    const banner = screen.getByRole("alert");
    expect(banner).toBeInTheDocument();
    expect(banner.textContent).toMatch(/could not build the series/i);
  });
});
