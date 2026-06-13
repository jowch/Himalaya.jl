import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, act, within } from "@testing-library/react";
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
// Records the exposure ids the page actually asks traces for (the fan-out) —
// lets the candidate-preview cap pin that below-fold candidates are not fetched.
let requestedTraceExposureIds: number[] = [];
// Per-exposure phase override for useMemberIndices (default "Pn3m") — lets
// preview-strip tests tell members apart by their accessible phase labels.
let phaseByExposure: Record<number, string> = {};

vi.mock("../../src/queries", () => ({
  useCorpusSampleTags: () => tagsState,
  useCorpusPickerSamples: () => pickerState,
  useScopeSeries: () => scopeState,
  useCreateSeries: () => createState,
  useMemberTraces: (ids: number[]) => {
    requestedTraceExposureIds = ids;
    return new Map<number, Trace>(ids.map((id) => [id, trace()]));
  },
  useMemberIndices: (ids: number[]) =>
    new Map<number, IndexEntry[]>(
      ids.map((id) => [id, [indexEntry(id, phaseByExposure[id] ?? "Pn3m", 0.9)]]),
    ),
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

/** A corpus exposing TWO distinct ordering variables on shared samples.
 * "ratio" is the frequency winner (3 pairs vs 2); "temp" is the alternative.
 * Each sample carries a distinct value under each key so a re-group is
 * observable in the rendered member values. */
function seedTwoKeys(): void {
  tagsState = {
    data: [
      { key: "ratio", value: "1 : 0" },
      { key: "ratio", value: "1 : 0.5" },
      { key: "ratio", value: "1 : 1" },
      { key: "temp", value: "20C" },
      { key: "temp", value: "40C" },
    ],
    isLoading: false,
    isError: false,
  };
  pickerState = {
    data: [
      pickerRow(sample(1, "A", [tag("ratio", "1 : 0"), tag("temp", "40C")]), 37),
      pickerRow(sample(2, "B", [tag("ratio", "1 : 0.5"), tag("temp", "20C")]), 65),
    ],
    isLoading: false,
    isError: false,
  };
}

/** 2 members + 7 loose candidates (no "ratio" tag) — loose exceeds the
 * candidate preview cap, so the worksheet must truncate to exemplars. Seven
 * (not six) so hidden (4) ≠ visible (3): a remainder computed from the SLICED
 * list instead of the full one is detectably wrong. */
function seedManyLoose(): void {
  tagsState = { data: [{ key: "ratio", value: "1 : 0.5" }], isLoading: false, isError: false };
  pickerState = {
    data: [
      pickerRow(sample(1, "A", [tag("ratio", "1 : 0")]), 37),
      pickerRow(sample(2, "B", [tag("ratio", "1 : 0.5")]), 65),
      pickerRow(sample(11, "L1", []), 100),
      pickerRow(sample(12, "L2", []), 101),
      pickerRow(sample(13, "L3", []), 102),
      pickerRow(sample(14, "L4", []), 103),
      pickerRow(sample(15, "L5", []), 104),
      pickerRow(sample(16, "L6", []), 105),
      pickerRow(sample(17, "L7", []), 106),
    ],
    isLoading: false,
    isError: false,
  };
}

beforeEach(() => {
  vi.clearAllMocks();
  scopeState = { mutate: scopeMutate, isSuccess: false, error: null, data: undefined };
  createState = { mutate: createMutate, isSuccess: false, error: null, data: undefined };
  requestedTraceExposureIds = [];
  phaseByExposure = {};
  seed();
});

describe("SeriesScopingPage", () => {
  it("renders a scope-sample-row per member", () => {
    renderPage();
    expect(screen.getAllByTestId("scope-sample-row")).toHaveLength(2);
  });

  it("names the skip gesture in the intro copy (SC-SKIPDISC)", () => {
    renderPage();
    // The skip toggle is otherwise invisible until hover; the intro must tell
    // a first-time user HOW, not just that skipping exists.
    expect(
      screen.getByText(/click a value to skip that sample/i),
    ).toBeInTheDocument();
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

  // ── SC-CANDID: candidate rows carry the member-row sample identity ──────────
  // Display names duplicate in the real corpus; without the mono smp_{id} the
  // page reads as a contradiction (the same name both in and out of the series).

  it("gives every visible candidate row its smp_{id} identity as rendered text", () => {
    seedManyLoose();
    renderPage();
    const candidates = screen.getAllByTestId("scope-candidate");
    expect(candidates).toHaveLength(3);
    const ids = candidates.map(
      (row) => within(row).getByText(/^smp_\d+$/).textContent,
    );
    expect(ids).toEqual(["smp_11", "smp_12", "smp_13"]);
  });

  it("distinguishes a member and a candidate sharing a display name by their smp ids", () => {
    // Member "1-2 Only" = smp_5; a DIFFERENT sample with the same display name
    // lacks the ratio key and splits out as a loose candidate (smp_42).
    pickerState = {
      data: [
        pickerRow(sample(5, "1-2 Only", [tag("ratio", "1 : 0")]), 37),
        pickerRow(sample(42, "1-2 Only", []), 100),
      ],
      isLoading: false,
      isError: false,
    };
    renderPage();
    const member = screen.getByTestId("scope-sample-row");
    expect(within(member).getByText("1-2 Only")).toBeInTheDocument();
    expect(within(member).getByText("smp_5")).toBeInTheDocument();
    const candidate = screen.getByTestId("scope-candidate");
    expect(within(candidate).getByText("1-2 Only")).toBeInTheDocument();
    expect(within(candidate).getByText("smp_42")).toBeInTheDocument();
  });

  it("names each loose candidate's phase as text, not by sparkline colour alone (SC-MISC / WCAG 1.4.1)", () => {
    // The candidate's only phase channel was the aria-hidden sparkline hue —
    // invisible to a colour-blind or screen-reader user. Its phase must also
    // read as TEXT in the row.
    phaseByExposure = { 100: "Ia3d" };
    pickerState = {
      data: [
        pickerRow(sample(5, "Member", [tag("ratio", "1 : 0")]), 37),
        pickerRow(sample(42, "Loose one", []), 100),
      ],
      isLoading: false,
      isError: false,
    };
    renderPage();
    const candidate = screen.getByTestId("scope-candidate");
    expect(within(candidate).getByText("Ia3d")).toBeInTheDocument();
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

  it("Discard sits in the plate foot beside Confirm & build, not at the page top (SC-POLISH2)", () => {
    seed3();
    renderPage();
    const discard = screen.getByTestId("scoping-discard");
    const build = screen.getByRole("button", { name: /confirm & build/i });
    // Same foot action cluster → they share a parent.
    expect(discard.parentElement).toBe(build.parentElement);
  });

  it("Discard on a clean worksheet leaves immediately with no confirm (SC-POLISH2)", () => {
    seed3();
    renderPage();
    fireEvent.click(screen.getByTestId("scoping-discard"));
    expect(navigateSpy).toHaveBeenCalledWith("/series");
    expect(screen.queryByTestId("scoping-discard-confirm")).not.toBeInTheDocument();
  });

  it("Discard after an edit asks first and never destroys the work until confirmed (SC-POLISH2)", () => {
    seed3();
    renderPage();
    // Dirty the worksheet: skip a member (pushes a history entry).
    fireEvent.click(screen.getAllByTestId("flag-button")[1]!);
    // Discard now arms an inline confirm — no navigation yet.
    fireEvent.click(screen.getByTestId("scoping-discard"));
    expect(screen.getByTestId("scoping-discard-confirm")).toBeInTheDocument();
    expect(navigateSpy).not.toHaveBeenCalledWith("/series");
    // "Keep editing" backs out, work intact.
    fireEvent.click(screen.getByRole("button", { name: /keep editing/i }));
    expect(screen.queryByTestId("scoping-discard-confirm")).not.toBeInTheDocument();
    expect(navigateSpy).not.toHaveBeenCalledWith("/series");
    // Re-arm and confirm → only now does it leave.
    fireEvent.click(screen.getByTestId("scoping-discard"));
    fireEvent.click(screen.getByTestId("scoping-discard-yes"));
    expect(navigateSpy).toHaveBeenCalledWith("/series");
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

  // ── SC-PREVIEWSKIP: the strip previews exactly the commit membership ───────
  // The PhaseStrip's job is to preview the series figure that will be built —
  // it must show exactly the members the write will commit (isKept), in
  // displayed order, never the skipped ones the foot says are excluded.

  it("the preview strip omits a skipped member's cell and restores it on unskip", () => {
    seed3();
    // Make the middle member (B → exposure 65) distinguishable by label.
    phaseByExposure = { 65: "Ia3d" };
    renderPage();
    const labels = (): (string | null)[] =>
      screen.getAllByTestId("ps-seg").map((el) => el.getAttribute("aria-label"));
    // All three members preview, in displayed (low→high) order.
    expect(labels()).toEqual(["Pn3m", "Ia3d", "Pn3m"]);
    // Skip B → its cell leaves the preview (2 cells, no Ia3d label).
    fireEvent.click(screen.getAllByTestId("flag-button")[1]!);
    expect(labels()).toEqual(["Pn3m", "Pn3m"]);
    // Unskip → the cell returns in place.
    fireEvent.click(screen.getAllByTestId("flag-button")[1]!);
    expect(labels()).toEqual(["Pn3m", "Ia3d", "Pn3m"]);
  });

  it("the preview strip follows the DISPLAYED order after a keyboard reorder, filtered by isKept", () => {
    // seed3's picker (rows) order and the value-sorted display order coincide,
    // so the omit-on-skip pin alone cannot tell `sorted` from `rows` — diverge
    // them with a real reorder (the SC-KBD grip path) before asserting.
    seed3();
    // Distinct label per member: A(37)=Pn3m (default), B(65)=Ia3d, C(66)=Im3m.
    phaseByExposure = { 65: "Ia3d", 66: "Im3m" };
    renderPage();
    const labels = (): (string | null)[] =>
      screen.getAllByTestId("ps-seg").map((el) => el.getAttribute("aria-label"));
    expect(labels()).toEqual(["Pn3m", "Ia3d", "Im3m"]);
    // ArrowUp on C's grip → displayed order A, C, B (rows order stays A, B, C).
    fireEvent.keyDown(screen.getByRole("button", { name: /^reorder C$/i }), { key: "ArrowUp" });
    expect(labels()).toEqual(["Pn3m", "Im3m", "Ia3d"]);
    // Skip A (first DISPLAYED row) → the kept set IN the new order: C, B.
    fireEvent.click(screen.getAllByTestId("flag-button")[0]!);
    expect(labels()).toEqual(["Im3m", "Ia3d"]);
  });

  it("omits the preview strip region entirely when every member is skipped", () => {
    renderPage();
    fireEvent.click(screen.getAllByTestId("flag-button")[0]!);
    fireEvent.click(screen.getAllByTestId("flag-button")[1]!);
    // Nothing will be committed → nothing to preview: no cells, no section
    // heading, no orphaned "No clear phase" caption. The foot already warns.
    expect(screen.queryAllByTestId("ps-seg")).toHaveLength(0);
    expect(
      screen.queryByRole("heading", { level: 2, name: /preview · phase across the series/i }),
    ).not.toBeInTheDocument();
    expect(screen.queryByText(/no clear phase/i)).not.toBeInTheDocument();
  });

  it("⌘Z from the page body steps the last skip back (regression pin)", () => {
    seed3();
    renderPage();
    fireEvent.click(screen.getAllByTestId("flag-button")[1]!);
    expect(screen.getByText(/2 values ready to commit · 1 skipped/i)).toBeInTheDocument();
    const notPrevented = fireEvent.keyDown(window, { key: "z", metaKey: true });
    // The page handled (and preventDefault'ed) the chord.
    expect(notPrevented).toBe(false);
    expect(screen.getByText(/3 values ready to commit/i)).toBeInTheDocument();
  });

  it("⌘Z from a text input does NOT fire skip-undo and leaves native text-undo intact", () => {
    seed3();
    renderPage();
    fireEvent.click(screen.getAllByTestId("flag-button")[1]!);
    expect(screen.getByText(/2 values ready to commit · 1 skipped/i)).toBeInTheDocument();
    const input = document.createElement("input");
    document.body.appendChild(input);
    try {
      const notPrevented = fireEvent.keyDown(input, { key: "z", metaKey: true });
      // Not preventDefault'ed — the browser's native text undo survives.
      expect(notPrevented).toBe(true);
      // The skip state is untouched.
      expect(screen.getByText(/2 values ready to commit · 1 skipped/i)).toBeInTheDocument();
    } finally {
      input.remove();
    }
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

  it("shows the tag-write failure copy when Op A fails (series was not created)", () => {
    scopeState = { mutate: scopeMutate, isSuccess: false, error: new Error("boom"), data: undefined };
    renderPage();
    const banner = screen.getByTestId("scoping-error-banner");
    expect(banner).toHaveAttribute("role", "alert");
    expect(banner.textContent).toMatch(/could not save the ordering tags/i);
    expect(banner.textContent).toMatch(/the series was not created/i);
    // The Op-B copy must NOT bleed in.
    expect(banner.textContent).not.toMatch(/were saved/i);
  });

  it("shows the tags-committed copy when Op B fails (Op A already wrote the tags)", () => {
    createState = { mutate: createMutate, isSuccess: false, error: new Error("boom"), data: undefined };
    renderPage();
    const banner = screen.getByTestId("scoping-error-banner");
    expect(banner).toHaveAttribute("role", "alert");
    expect(banner.textContent).toMatch(/the ordering tags were saved, but the series was not created/i);
    // Op A committed — the page must never claim nothing was saved.
    expect(screen.queryByText(/nothing was saved/i)).not.toBeInTheDocument();
    // The Op-A copy must NOT bleed in.
    expect(banner.textContent).not.toMatch(/could not save/i);
  });

  it("prefers the Op-A copy when both errors are set (a stale create error can linger from a prior attempt)", () => {
    scopeState = { mutate: scopeMutate, isSuccess: false, error: new Error("fresh"), data: undefined };
    createState = { mutate: createMutate, isSuccess: false, error: new Error("stale"), data: undefined };
    renderPage();
    const banner = screen.getByTestId("scoping-error-banner");
    expect(banner.textContent).toMatch(/could not save the ordering tags/i);
    expect(banner.textContent).not.toMatch(/were saved/i);
  });

  it("SC-POLISH2: Confirm & build reads Building… with aria-busy while the scope→create chain runs", () => {
    const { rerender } = renderPage();
    fireEvent.click(screen.getByRole("button", { name: /confirm & build/i }));
    // Op A is in flight → the control tells the truth: progressive register +
    // aria-busy, still disabled (no double-submit).
    act(() => rerender());
    const busy = screen.getByRole("button", { name: /building…/i });
    expect(busy).toHaveAttribute("aria-busy", "true");
    expect(busy).toBeDisabled();
    expect(screen.queryByRole("button", { name: /confirm & build/i })).not.toBeInTheDocument();

    // Op A lands; the gap before Op B settles is part of the chain → still busy.
    scopeState = { ...scopeState, isSuccess: true };
    act(() => rerender());
    expect(createMutate).toHaveBeenCalled();
    expect(screen.getByRole("button", { name: /building…/i })).toHaveAttribute(
      "aria-busy",
      "true",
    );
  });

  it("SC-POLISH2: the busy register reverts when Op A fails (truthful-error exit, retry re-armed)", () => {
    const { rerender } = renderPage();
    fireEvent.click(screen.getByRole("button", { name: /confirm & build/i }));
    act(() => rerender());
    expect(screen.getByRole("button", { name: /building…/i })).toBeInTheDocument();

    // The tag write FAILS. The render that surfaces the error banner must
    // already show the resting register — the stage ref resets in an effect
    // that triggers no re-render of its own.
    scopeState = { ...scopeState, error: new Error("boom") };
    act(() => rerender());
    expect(screen.getByTestId("scoping-error-banner")).toBeInTheDocument();
    const resting = screen.getByRole("button", { name: /confirm & build/i });
    expect(resting).not.toHaveAttribute("aria-busy");
    expect(resting).not.toBeDisabled();
    expect(screen.queryByRole("button", { name: /building…/i })).not.toBeInTheDocument();
  });

  it("SC-POLISH2: the busy register reverts when Op B fails (tags written, create missing)", () => {
    const { rerender } = renderPage();
    fireEvent.click(screen.getByRole("button", { name: /confirm & build/i }));
    scopeState = { ...scopeState, isSuccess: true };
    act(() => rerender());
    expect(createMutate).toHaveBeenCalled();
    expect(screen.getByRole("button", { name: /building…/i })).toBeInTheDocument();

    createState = { ...createState, error: new Error("boom") };
    act(() => rerender());
    expect(screen.getByTestId("scoping-error-banner").textContent).toMatch(
      /the ordering tags were saved/i,
    );
    const resting = screen.getByRole("button", { name: /confirm & build/i });
    expect(resting).not.toHaveAttribute("aria-busy");
    expect(screen.queryByRole("button", { name: /building…/i })).not.toBeInTheDocument();
  });

  it("SC-POLISH2: the cold-path foot button carries the same busy register through the same chain", () => {
    tagsState = { data: [], isLoading: false, isError: false };
    pickerState = {
      data: [pickerRow(sample(1, "A", []), 37), pickerRow(sample(2, "B", []), 65)],
      isLoading: false,
      isError: false,
    };
    const { rerender } = renderPage([1, 2]);
    fireEvent.change(screen.getByTestId("cold-key-input"), { target: { value: "dose" } });
    const valueInputs = screen.getAllByTestId("cold-assign-row");
    fireEvent.change(valueInputs[0]!.querySelector("input")!, { target: { value: "10" } });
    fireEvent.change(valueInputs[1]!.querySelector("input")!, { target: { value: "20" } });
    fireEvent.click(screen.getByRole("button", { name: /confirm & build/i }));

    act(() => rerender());
    const busy = screen.getByRole("button", { name: /building…/i });
    expect(busy).toHaveAttribute("aria-busy", "true");
    expect(busy).toBeDisabled();

    // Op A fails → reverts to the resting register, re-armed for a retry.
    scopeState = { ...scopeState, error: new Error("boom") };
    act(() => rerender());
    const resting = screen.getByRole("button", { name: /confirm & build/i });
    expect(resting).not.toHaveAttribute("aria-busy");
    expect(resting).not.toBeDisabled();
  });

  it("renders the ordered-by control as a dropdown when ≥2 ordering variables exist", () => {
    seedTwoKeys();
    renderPage();
    const trigger = screen.getByTestId("order-field");
    // Dropdown mode → the trigger advertises a menu popup (static mode has none).
    expect(trigger).toHaveAttribute("aria-haspopup", "menu");
  });

  it("lists a menu option for each distinct ordering variable", () => {
    seedTwoKeys();
    renderPage();
    fireEvent.click(screen.getByTestId("order-field"));
    expect(screen.getByRole("menuitemradio", { name: "ratio" })).toBeInTheDocument();
    expect(screen.getByRole("menuitemradio", { name: "temp" })).toBeInTheDocument();
    // Both corpus keys plus the always-present "Define your own…" sentinel.
    expect(screen.getAllByRole("menuitemradio")).toHaveLength(3);
  });

  it("lists Define your own… as the last dropdown option", () => {
    seedTwoKeys();
    renderPage();
    fireEvent.click(screen.getByTestId("order-field"));
    const items = screen.getAllByRole("menuitemradio").map((m) => m.textContent);
    expect(items[items.length - 1]).toBe("Define your own…");
  });

  it("selecting a non-active ordering variable re-groups the worksheet by its values", () => {
    seedTwoKeys();
    renderPage();
    // Default: grouped by the frequency winner "ratio" → ratio values render.
    // (FlagButton textContent carries an sr-only "Skip this read: " purpose prefix.)
    let values = screen.getAllByTestId("flag-button").map((b) => b.textContent ?? "");
    expect(values).toEqual(
      expect.arrayContaining(["Skip this read: 1 : 0", "Skip this read: 1 : 0.5"]),
    );
    expect(values.some((v) => v.includes("20C"))).toBe(false);

    // Switch the ordering variable to "temp".
    fireEvent.click(screen.getByTestId("order-field"));
    fireEvent.click(screen.getByRole("menuitemradio", { name: "temp" }));

    // The members now read their "temp" values.
    values = screen.getAllByTestId("flag-button").map((b) => b.textContent ?? "");
    expect(values).toEqual(
      expect.arrayContaining(["Skip this read: 20C", "Skip this read: 40C"]),
    );
    expect(values.some((v) => v.includes("1 : 0"))).toBe(false);
  });

  it("renders the dropdown even when only one ordering variable exists (Define your own… is always an alternative)", () => {
    // Deliberately flipped from the commit-1 static-control behaviour: with the
    // "Define your own…" sentinel there is ALWAYS a real alternative, so the
    // single-key `seed()` corpus now gets the dropdown too.
    renderPage();
    const field = screen.getByTestId("order-field");
    expect(field).toHaveAttribute("aria-haspopup", "menu");
    fireEvent.click(field);
    const items = screen.getAllByRole("menuitemradio").map((m) => m.textContent);
    expect(items).toEqual(["ratio", "Define your own…"]);
  });

  it("selecting Define your own… swaps the worksheet for the custom assign card seeded with the members", () => {
    seedTwoKeys();
    renderPage();
    fireEvent.click(screen.getByTestId("order-field"));
    fireEvent.click(screen.getByRole("menuitemradio", { name: "Define your own…" }));

    // The warm worksheet is replaced by the custom card in the same slot.
    expect(screen.getByTestId("custom-scope-plate")).toBeInTheDocument();
    expect(screen.queryByTestId("scope-sample-row")).toBeNull();
    expect(screen.getByTestId("cold-assign-panel")).toBeInTheDocument();

    // Headed, not headless: the h1 reads a placeholder until the key is named.
    expect(screen.getByRole("heading", { level: 1 }).textContent).toBe("Series by …");

    // One assign row per current member, names carried over, values empty.
    const rows = screen.getAllByTestId("cold-assign-row");
    expect(rows).toHaveLength(2);
    expect(rows[0]!.textContent).toContain("A");
    expect(rows[1]!.textContent).toContain("B");
    for (const row of rows) {
      expect((row.querySelector("input") as HTMLInputElement).value).toBe("");
    }
  });

  it("custom mode suppresses the cold-corpus intro copy (these samples DO share keys)", () => {
    seedTwoKeys();
    renderPage();
    fireEvent.click(screen.getByTestId("order-field"));
    fireEvent.click(screen.getByRole("menuitemradio", { name: "Define your own…" }));
    expect(screen.getByTestId("custom-scope-plate")).toBeInTheDocument();
    // The default cold-corpus paragraph would be false here and is redundant
    // with the card's own caption, so it must not render.
    expect(screen.queryByText(/these samples share no tag key yet/i)).toBeNull();
    // The card's own instruction stands in for it.
    expect(screen.getByText(/name the variable below and assign each sample's value/i)).toBeInTheDocument();
  });

  it("the cold path keeps the cold-corpus intro copy (it is true there)", () => {
    // No proposable key + a deliberate seed → the cold path, where the samples
    // genuinely share no tag key, so the default intro must still render.
    tagsState = { data: [], isLoading: false, isError: false };
    pickerState = {
      data: [pickerRow(sample(1, "A", []), 37), pickerRow(sample(2, "B", []), 65)],
      isLoading: false,
      isError: false,
    };
    renderPage([1, 2]);
    expect(screen.getByTestId("cold-scope-plate")).toBeInTheDocument();
    expect(screen.getByText(/these samples share no tag key yet/i)).toBeInTheDocument();
  });

  it("custom mode with zero members keeps Confirm & build disabled even after typing a key", () => {
    // A key IS proposed from the corpus-wide tag pairs, but no picker sample
    // carries it → a warm worksheet with zero members (all loose candidates).
    // Entering custom mode then seeds ZERO assign rows; a key alone must not
    // open the gate (an empty worksheet never builds an empty series).
    pickerState = {
      data: [pickerRow(sample(1, "A", []), 37), pickerRow(sample(2, "B", []), 65)],
      isLoading: false,
      isError: false,
    };
    renderPage();
    fireEvent.click(screen.getByTestId("order-field"));
    fireEvent.click(screen.getByRole("menuitemradio", { name: "Define your own…" }));
    expect(screen.getByTestId("custom-scope-plate")).toBeInTheDocument();
    expect(screen.queryAllByTestId("cold-assign-row")).toHaveLength(0);
    fireEvent.change(screen.getByTestId("cold-key-input"), { target: { value: "dose" } });
    expect(screen.getByRole("button", { name: /confirm & build/i })).toBeDisabled();
  });

  it("custom mode: types a key, assigns values, runs the same scope-then-create chain, navigates to /series/:id", () => {
    seedTwoKeys();
    const { rerender } = renderPage();
    fireEvent.click(screen.getByTestId("order-field"));
    fireEvent.click(screen.getByRole("menuitemradio", { name: "Define your own…" }));

    // Gate closed until the key and every value are filled.
    expect(screen.getByRole("button", { name: /confirm & build/i })).toBeDisabled();
    fireEvent.change(screen.getByTestId("cold-key-input"), { target: { value: "dose" } });
    const rows = screen.getAllByTestId("cold-assign-row");
    fireEvent.change(rows[0]!.querySelector("input")!, { target: { value: "10" } });
    fireEvent.change(rows[1]!.querySelector("input")!, { target: { value: "20" } });
    // The heading follows the typed key.
    expect(screen.getByRole("heading", { level: 1 }).textContent).toBe("Series by dose");
    const build = screen.getByRole("button", { name: /confirm & build/i });
    expect(build).not.toBeDisabled();
    fireEvent.click(build);

    // Op A — the tag batch write under the typed key, values in worksheet order.
    expect(scopeMutate).toHaveBeenCalledWith(
      expect.objectContaining({
        key: "dose",
        tags: [
          { sampleId: 1, value: "10" },
          { sampleId: 2, value: "20" },
        ],
      }),
    );
    expect(createMutate).not.toHaveBeenCalled();

    // Op A lands → Op B fires with the create body built from the custom rows.
    scopeState = { mutate: scopeMutate, isSuccess: true, error: null, data: undefined };
    act(() => rerender());
    expect(createMutate).toHaveBeenCalledWith(
      expect.objectContaining({
        title: "Series by dose",
        ordering_variable: "dose",
        samples: [
          { sample_id: 1, position: 0 },
          { sample_id: 2, position: 1 },
        ],
      }),
    );

    // Op B lands with the created Series → navigate to the new builder.
    createState = { mutate: createMutate, isSuccess: true, error: null, data: fullSeries(11) };
    act(() => rerender());
    expect(navigateSpy).toHaveBeenCalledWith("/series/11");
  });

  // ── SC-FOLD: loose-candidate preview cap ────────────────────────────────────
  // The worksheet is a review surface; on a whole-corpus visit the loose list
  // can run to 100+ rows, burying the members and the build action. The page
  // renders only a few exemplar candidates plus one honest section note.

  it("caps the rendered candidate rows at the preview count", () => {
    seedManyLoose();
    renderPage();
    // 7 loose candidates, but only the preview cap (3) render as rows.
    expect(screen.getAllByTestId("scope-candidate")).toHaveLength(3);
    // Members are unaffected.
    expect(screen.getAllByTestId("scope-sample-row")).toHaveLength(2);
  });

  it("the section note states the hidden remainder count", () => {
    seedManyLoose();
    renderPage();
    const note = screen.getByTestId("scope-candidates-note");
    // 7 loose − 3 shown = 4 more (computed from the FULL list; 4 ≠ the visible
    // count, so a remainder derived from the sliced list fails here).
    expect(note.textContent).toMatch(/4 more lack the ratio/i);
  });

  it("pluralizes a single hidden candidate as 'lacks'", () => {
    // 1 member + 4 loose: exactly one is hidden behind the cap of 3.
    pickerState = {
      data: [
        pickerRow(sample(1, "A", [tag("ratio", "1 : 0")]), 37),
        pickerRow(sample(11, "L1", []), 100),
        pickerRow(sample(12, "L2", []), 101),
        pickerRow(sample(13, "L3", []), 102),
        pickerRow(sample(14, "L4", []), 103),
      ],
      isLoading: false,
      isError: false,
    };
    renderPage();
    const note = screen.getByTestId("scope-candidates-note");
    expect(note.textContent).toMatch(/1 more lacks the ratio/i);
  });

  it("the note's contact-sheet control navigates to /samples", () => {
    seedManyLoose();
    renderPage();
    const note = screen.getByTestId("scope-candidates-note");
    fireEvent.click(within(note).getByRole("button", { name: /open the contact sheet/i }));
    expect(navigateSpy).toHaveBeenCalledWith("/samples");
  });

  it("renders every candidate and no remainder line when loose fits within the cap", () => {
    // 1 member + 2 loose (under the cap of 3): nothing is hidden, so the note
    // must not claim a remainder — but the instruction + control still render.
    pickerState = {
      data: [
        pickerRow(sample(1, "A", [tag("ratio", "1 : 0")]), 37),
        pickerRow(sample(11, "L1", []), 100),
        pickerRow(sample(12, "L2", []), 101),
      ],
      isLoading: false,
      isError: false,
    };
    renderPage();
    expect(screen.getAllByTestId("scope-candidate")).toHaveLength(2);
    const note = screen.getByTestId("scope-candidates-note");
    expect(note.textContent).not.toMatch(/more lack/i);
    expect(
      within(note).getByRole("button", { name: /open the contact sheet/i }),
    ).toBeInTheDocument();
  });

  it("keeps the nothing-else-matches line when there are no candidates", () => {
    // seed(): both samples carry the key → loose is empty. The existing branch
    // stands, and the candidates note must not render.
    renderPage();
    expect(
      screen.getByText(/nothing else in the corpus matches this grouping/i),
    ).toBeInTheDocument();
    expect(screen.queryByTestId("scope-candidates-note")).toBeNull();
  });

  it("fetches traces only for members plus the visible candidates", () => {
    seedManyLoose();
    renderPage();
    // 2 member exposures (37, 65) + the 3 visible candidates' (100–102).
    expect(requestedTraceExposureIds).toHaveLength(5);
    expect(requestedTraceExposureIds).toEqual(expect.arrayContaining([37, 65]));
    // The hidden candidates' exposures must NOT be fan-out fetched.
    expect(requestedTraceExposureIds).not.toContain(103);
    expect(requestedTraceExposureIds).not.toContain(104);
    expect(requestedTraceExposureIds).not.toContain(105);
    expect(requestedTraceExposureIds).not.toContain(106);
  });

  it("selecting an existing key from custom mode returns to the warm worksheet re-grouped by it", () => {
    seedTwoKeys();
    renderPage();
    fireEvent.click(screen.getByTestId("order-field"));
    fireEvent.click(screen.getByRole("menuitemradio", { name: "Define your own…" }));
    expect(screen.getByTestId("custom-scope-plate")).toBeInTheDocument();

    // The custom card's order field is the escape hatch back to a proposed key.
    fireEvent.click(screen.getByTestId("order-field"));
    fireEvent.click(screen.getByRole("menuitemradio", { name: "temp" }));

    expect(screen.queryByTestId("custom-scope-plate")).toBeNull();
    const values = screen.getAllByTestId("flag-button").map((b) => b.textContent ?? "");
    expect(values).toEqual(
      expect.arrayContaining(["Skip this read: 20C", "Skip this read: 40C"]),
    );
    expect(values.some((v) => v.includes("1 : 0"))).toBe(false);
  });

  // ── SC-SEEDDEAD: the proposal vocabulary scopes to the seed ────────────────
  // A contact-sheet seed of samples lacking the corpus's dominant tag key must
  // never dead-end: the ordering-variable universe is the SEEDED samples' own
  // tags, not the corpus-wide pair list (which applies only to direct visits).
  describe("seed-scoped vocabulary (SC-SEEDDEAD)", () => {
    /** Corpus where "ratio" is the frequency winner; samples 11-13 untagged. */
    function seedUntaggedCorpus(): void {
      tagsState = {
        data: [
          { key: "ratio", value: "1 : 0" },
          { key: "ratio", value: "1 : 0.5" },
        ],
        isLoading: false,
        isError: false,
      };
      pickerState = {
        data: [
          pickerRow(sample(1, "A", [tag("ratio", "1 : 0")]), 37),
          pickerRow(sample(2, "B", [tag("ratio", "1 : 0.5")]), 65),
          pickerRow(sample(11, "U1", []), 100),
          pickerRow(sample(12, "U2", []), 101),
          pickerRow(sample(13, "U3", []), 102),
        ],
        isLoading: false,
        isError: false,
      };
    }

    /** Corpus winner is "ratio" (3 pairs vs 2); samples 4-5 carry ONLY "temp". */
    function seedDisjointKeyCorpus(): void {
      tagsState = {
        data: [
          { key: "ratio", value: "1 : 0" },
          { key: "ratio", value: "1 : 0.5" },
          { key: "ratio", value: "1 : 1" },
          { key: "temp", value: "20C" },
          { key: "temp", value: "40C" },
        ],
        isLoading: false,
        isError: false,
      };
      pickerState = {
        data: [
          pickerRow(sample(1, "A", [tag("ratio", "1 : 0")]), 37),
          pickerRow(sample(2, "B", [tag("ratio", "1 : 0.5")]), 65),
          pickerRow(sample(3, "C", [tag("ratio", "1 : 1")]), 66),
          pickerRow(sample(4, "T1", [tag("temp", "20C")]), 70),
          pickerRow(sample(5, "T2", [tag("temp", "40C")]), 71),
        ],
        isLoading: false,
        isError: false,
      };
    }

    it("seeding only untagged samples on a tagged corpus renders the cold-assign card, not a 0-member worksheet", () => {
      // The corpus proposes "ratio", but no seeded sample carries ANY tag —
      // there is no variable proposable FROM THE SEED, so the cold-assign
      // escape must fire with every seeded sample as an assign row.
      seedUntaggedCorpus();
      renderPage([11, 12, 13]);
      expect(screen.getByTestId("cold-scope-plate")).toBeInTheDocument();
      const rows = screen.getAllByTestId("cold-assign-row");
      expect(rows).toHaveLength(3);
      expect(rows[0]!.textContent).toContain("U1");
      expect(rows[1]!.textContent).toContain("U2");
      expect(rows[2]!.textContent).toContain("U3");
      // The dead-end 0-member warm worksheet must NOT render.
      expect(screen.queryByTestId("scope-sample-row")).toBeNull();
      expect(screen.queryByText(/keep at least one value to build/i)).toBeNull();
    });

    it("proposes the seeded samples' own shared key, not the corpus-wide frequency winner", () => {
      // Corpus-wide "ratio" outweighs "temp", but the seeded samples carry
      // only "temp" → the worksheet must group by temp with both as members.
      seedDisjointKeyCorpus();
      renderPage([4, 5]);
      expect(screen.getByRole("heading", { level: 1 }).textContent).toBe("Series by temp");
      expect(screen.getAllByTestId("scope-sample-row")).toHaveLength(2);
    });

    it("offers only keys the seeded samples carry in the order dropdown (plus Define your own…)", () => {
      // A corpus key absent from the seed ("ratio") would re-create the
      // 0-member state straight from the dropdown — it must not be offered.
      seedDisjointKeyCorpus();
      renderPage([4, 5]);
      fireEvent.click(screen.getByTestId("order-field"));
      const items = screen.getAllByRole("menuitemradio").map((m) => m.textContent);
      expect(items).toEqual(["temp", "Define your own…"]);
    });

    it("Define your own… on a seeded visit carries the WHOLE selection (members first, then loose)", () => {
      // Mixed seed: 2 with "ratio" + 2 without. The warm worksheet shows
      // 2 members + 2 candidates; custom mode must seed all 4 (the user
      // deliberately selected them), members in displayed order first.
      seedUntaggedCorpus();
      renderPage([1, 2, 11, 12]);
      expect(screen.getAllByTestId("scope-sample-row")).toHaveLength(2);
      expect(screen.getAllByTestId("scope-candidate")).toHaveLength(2);
      fireEvent.click(screen.getByTestId("order-field"));
      fireEvent.click(screen.getByRole("menuitemradio", { name: "Define your own…" }));
      const rows = screen.getAllByTestId("cold-assign-row");
      expect(rows).toHaveLength(4);
      expect(rows[0]!.textContent).toContain("A");
      expect(rows[1]!.textContent).toContain("B");
      expect(rows[2]!.textContent).toContain("U1");
      expect(rows[3]!.textContent).toContain("U2");
    });

    it("Define your own… on a direct visit stays members-only (loose candidates do not fan in)", () => {
      // Whole-corpus DYO must NOT suck in the 7 loose candidates — only the
      // 2 members seed the custom worksheet.
      seedManyLoose();
      renderPage();
      fireEvent.click(screen.getByTestId("order-field"));
      fireEvent.click(screen.getByRole("menuitemradio", { name: "Define your own…" }));
      expect(screen.getAllByTestId("cold-assign-row")).toHaveLength(2);
    });
  });

  // ── SC-COLD + SC-COLDHEAD: cold-path plate identity + heading structure ────
  // The cold worksheet previously rendered a bare Card: no page h1 and its
  // "Name the ordering variable" label was a div Kicker (WCAG 1.3.1). It must
  // carry the SAME plate identity as the warm/custom plates ("New series"
  // kicker + an h1 that follows the typed key) and expose its section label as
  // a real h2 so the tree has no level skip.
  describe("cold-path identity and heading tree (SC-COLD / SC-COLDHEAD)", () => {
    /** No proposable key + a deliberate seed → the cold path. */
    function seedColdPath(): void {
      tagsState = { data: [], isLoading: false, isError: false };
      pickerState = {
        data: [pickerRow(sample(1, "A", []), 37), pickerRow(sample(2, "B", []), 65)],
        isLoading: false,
        isError: false,
      };
    }

    /** Asserts the rendered heading tree starts at h1 and never skips a level. */
    function expectNoHeadingLevelSkip(): void {
      const levels = screen.getAllByRole("heading").map((h) => Number(h.tagName.slice(1)));
      expect(Math.min(...levels)).toBe(1);
      for (const level of levels) {
        if (level > 1) expect(levels).toContain(level - 1);
      }
    }

    it("cold path renders the plate identity: New series kicker + an h1 that follows the typed key", () => {
      seedColdPath();
      renderPage([1, 2]);
      expect(screen.getByTestId("cold-scope-plate")).toBeInTheDocument();
      expect(screen.getByText("New series")).toBeInTheDocument();
      // Placeholder until the key is named (mirrors the custom card's h1).
      expect(screen.getByRole("heading", { level: 1 }).textContent).toBe("Series by …");
      fireEvent.change(screen.getByTestId("cold-key-input"), {
        target: { value: "lipid ratio" },
      });
      expect(screen.getByRole("heading", { level: 1 }).textContent).toBe(
        "Series by lipid ratio",
      );
    });

    it("cold path exposes the section label as an h2 under a single h1, with no level skip", () => {
      seedColdPath();
      renderPage([1, 2]);
      expect(screen.getAllByRole("heading", { level: 1 })).toHaveLength(1);
      expect(
        screen.getByRole("heading", { level: 2, name: "Name the ordering variable" }),
      ).toBeInTheDocument();
      expectNoHeadingLevelSkip();
    });

    it("custom (Define your own…) mode keeps its identity block and slots the same h2 section label in, no level skip", () => {
      seedTwoKeys();
      renderPage();
      fireEvent.click(screen.getByTestId("order-field"));
      fireEvent.click(screen.getByRole("menuitemradio", { name: "Define your own…" }));
      // The identity block is unchanged (pinned elsewhere too: placeholder h1 +
      // suppressed cold-corpus intro); here the level structure is the pin.
      expect(screen.getByRole("heading", { level: 1 }).textContent).toBe("Series by …");
      expect(screen.getByRole("heading", { level: 2, name: "Ordered by" })).toBeInTheDocument();
      expect(
        screen.getByRole("heading", { level: 2, name: "Name the ordering variable" }),
      ).toBeInTheDocument();
      expectNoHeadingLevelSkip();
    });

    it("the warm worksheet tree is unchanged: one h1, h2 section labels, no cold-assign heading", () => {
      seed3();
      renderPage();
      expect(screen.getAllByRole("heading", { level: 1 })).toHaveLength(1);
      expect(screen.getByRole("heading", { level: 2, name: "Ordered by" })).toBeInTheDocument();
      expect(
        screen.queryByRole("heading", { name: "Name the ordering variable" }),
      ).not.toBeInTheDocument();
      expectNoHeadingLevelSkip();
    });
  });

  describe("keyboard reorder (SC-KBD)", () => {
    /** Member names (A/B/C) in rendered worksheet order. */
    function renderedNames(): string[] {
      return screen
        .getAllByTestId("scope-sample-row")
        .map((r) => within(r).getByText(/^(A|B|C)$/).textContent ?? "");
    }

    it("ArrowUp on a middle row's grip button moves it up in the rendered order", () => {
      seed3();
      renderPage();
      expect(renderedNames()).toEqual(["A", "B", "C"]);
      fireEvent.keyDown(screen.getByRole("button", { name: /^reorder B$/i }), { key: "ArrowUp" });
      expect(renderedNames()).toEqual(["B", "A", "C"]);
    });

    it("announces a successful move with the computed position and count", () => {
      seed3();
      renderPage();
      fireEvent.keyDown(screen.getByRole("button", { name: /^reorder B$/i }), { key: "ArrowUp" });
      expect(screen.getByTestId("reorder-announcement")).toHaveTextContent(
        "Moved B to position 1 of 3.",
      );
    });

    it("ArrowUp on the first row is an order no-op but still announces 'already first'", () => {
      seed3();
      renderPage();
      fireEvent.keyDown(screen.getByRole("button", { name: /^reorder A$/i }), { key: "ArrowUp" });
      expect(renderedNames()).toEqual(["A", "B", "C"]);
      expect(screen.getByTestId("reorder-announcement")).toHaveTextContent("A is already first.");
    });

    it("ArrowDown on the last row is an order no-op but still announces 'already last'", () => {
      seed3();
      renderPage();
      fireEvent.keyDown(screen.getByRole("button", { name: /^reorder C$/i }), { key: "ArrowDown" });
      expect(renderedNames()).toEqual(["A", "B", "C"]);
      expect(screen.getByTestId("reorder-announcement")).toHaveTextContent("C is already last.");
    });

    it("keeps focus on the moved row's grip button after the move", () => {
      seed3();
      renderPage();
      const grip = screen.getByRole("button", { name: /^reorder B$/i });
      grip.focus();
      fireEvent.keyDown(grip, { key: "ArrowDown" });
      expect(renderedNames()).toEqual(["A", "C", "B"]);
      expect(document.activeElement).toBe(grip);
      // The focused node must still BE row B's grip (catches index-keyed rows,
      // where focus would silently land on a different sample's button).
      expect(screen.getByRole("button", { name: /^reorder B$/i })).toBe(document.activeElement);
      // A second computed announcement (resists a hardcoded string).
      expect(screen.getByTestId("reorder-announcement")).toHaveTextContent(
        "Moved B to position 3 of 3.",
      );
    });

    it("renders the polite live region persistently, before any move", () => {
      seed3();
      renderPage();
      const region = screen.getByTestId("reorder-announcement");
      expect(region).toHaveAttribute("aria-live", "polite");
      expect(region).toHaveTextContent("");
    });

    it("two identical consecutive announcements still mutate textContent (flip)", () => {
      // SRs skip a live-region update whose text is byte-identical to the last;
      // the toggled trailing space (LiveRegion.tsx pattern) must survive. Raw
      // textContent comparison on purpose: toHaveTextContent normalizes
      // whitespace and cannot see the flip.
      seed3();
      renderPage();
      const grip = screen.getByRole("button", { name: /^reorder A$/i });
      fireEvent.keyDown(grip, { key: "ArrowUp" });
      const first = screen.getByTestId("reorder-announcement").textContent;
      fireEvent.keyDown(grip, { key: "ArrowUp" });
      const second = screen.getByTestId("reorder-announcement").textContent;
      expect(first).toMatch(/A is already first\./);
      expect(second).toMatch(/A is already first\./);
      expect(second).not.toBe(first);
    });
  });

  // SC-TAGORDER: Op A (the tag write) must consume the SAME displayed/reordered
  // member order as Op B (the series position-create), so the recipe order and
  // the tag-write order never diverge. Pre-fix the tags came from the
  // never-reordered proposal `rows`; the fix routes both through `sorted`.
  describe("tag-write order follows the displayed order (SC-TAGORDER)", () => {
    it("writes the ordering tags in the reordered display order, matching the series positions", () => {
      seed3();
      const { rerender } = renderPage();
      // Default low→high display order is A, B, C (ids 1, 2, 3).
      // ArrowUp on C's grip → displayed order A, C, B (ids 1, 3, 2).
      fireEvent.keyDown(screen.getByRole("button", { name: /^reorder C$/i }), { key: "ArrowUp" });
      fireEvent.click(screen.getByRole("button", { name: /confirm & build/i }));
      // Op A tags follow the DISPLAYED order (exact array), not the proposal order.
      expect(scopeMutate).toHaveBeenCalledWith(
        expect.objectContaining({
          key: "ratio",
          tags: [
            { sampleId: 1, value: "1 : 0" },
            { sampleId: 3, value: "1 : 1" },
            { sampleId: 2, value: "1 : 0.5" },
          ],
        }),
      );
      // Op B positions follow the same order — write order and recipe agree.
      scopeState = { mutate: scopeMutate, isSuccess: true, error: null, data: undefined };
      act(() => rerender());
      expect(createMutate).toHaveBeenCalledWith(
        expect.objectContaining({
          samples: [
            { sample_id: 1, position: 0 },
            { sample_id: 3, position: 1 },
            { sample_id: 2, position: 2 },
          ],
        }),
      );
    });
  });

  // SC-COUNTHONEST: the count caption must not keep claiming "low to high" once
  // the order is no longer the value-sorted default (the controls-don't-lie law).
  describe("count caption honesty after reorder", () => {
    it("reads 'low to high' while the order is the canonical value-sorted default", () => {
      seed3();
      renderPage();
      expect(screen.getByText(/3 samples · low to high/i)).toBeInTheDocument();
      expect(screen.queryByText(/custom order/i)).toBeNull();
    });

    it("flips to 'custom order' after a manual reorder diverges from the default", () => {
      seed3();
      renderPage();
      fireEvent.keyDown(screen.getByRole("button", { name: /^reorder B$/i }), { key: "ArrowUp" });
      expect(screen.getByText(/3 samples · custom order/i)).toBeInTheDocument();
      expect(screen.queryByText(/low to high/i)).toBeNull();
    });

    it("returns to 'low to high' when the order is restored to the default", () => {
      seed3();
      renderPage();
      // Move B up, then back down — the order is the default again.
      fireEvent.keyDown(screen.getByRole("button", { name: /^reorder B$/i }), { key: "ArrowUp" });
      expect(screen.getByText(/custom order/i)).toBeInTheDocument();
      fireEvent.keyDown(screen.getByRole("button", { name: /^reorder B$/i }), { key: "ArrowDown" });
      expect(screen.getByText(/3 samples · low to high/i)).toBeInTheDocument();
    });
  });

  // SC-REORDERUNDO: a manual reorder is recoverable through the SAME undo path
  // (⌘Z + the plate's "Undo last change") that a skip is — no asymmetric freedom.
  describe("reorder undo (parity with skip)", () => {
    function renderedNames(): string[] {
      return screen
        .getAllByTestId("scope-sample-row")
        .map((r) => within(r).getByText(/^(A|B|C)$/).textContent ?? "");
    }

    it("⌘Z restores the prior display order after a reorder", () => {
      seed3();
      renderPage();
      expect(renderedNames()).toEqual(["A", "B", "C"]);
      fireEvent.keyDown(screen.getByRole("button", { name: /^reorder B$/i }), { key: "ArrowUp" });
      expect(renderedNames()).toEqual(["B", "A", "C"]);
      const notPrevented = fireEvent.keyDown(window, { key: "z", metaKey: true });
      expect(notPrevented).toBe(false);
      expect(renderedNames()).toEqual(["A", "B", "C"]);
      // And the caption is honest again.
      expect(screen.getByText(/3 samples · low to high/i)).toBeInTheDocument();
    });

    it("the plate's Undo affordance appears after a reorder and is labelled for it", () => {
      seed3();
      renderPage();
      // No edit yet → no Undo.
      expect(screen.queryByRole("button", { name: /undo/i })).toBeNull();
      fireEvent.keyDown(screen.getByRole("button", { name: /^reorder B$/i }), { key: "ArrowUp" });
      const undoBtn = screen.getByRole("button", { name: /undo last change/i });
      // Honest label for the reorder case (title attr; see ScopePlate undoLabel).
      expect(undoBtn.getAttribute("title")).toBe("Step back: reorder");
      // Clicking it restores the default order.
      fireEvent.click(undoBtn);
      expect(renderedNames()).toEqual(["A", "B", "C"]);
    });

    it("a boundary (no-op) reorder records no undo entry", () => {
      seed3();
      renderPage();
      // ArrowUp on the first row cannot move it — no history pushed, no Undo.
      fireEvent.keyDown(screen.getByRole("button", { name: /^reorder A$/i }), { key: "ArrowUp" });
      expect(screen.queryByRole("button", { name: /undo/i })).toBeNull();
      expect(screen.getByText(/3 samples · low to high/i)).toBeInTheDocument();
    });

    it("undo unwinds reorder and skip independently, last-in-first-out", () => {
      seed3();
      renderPage();
      // 1) reorder B up, 2) skip A (now at position 2).
      fireEvent.keyDown(screen.getByRole("button", { name: /^reorder B$/i }), { key: "ArrowUp" });
      expect(renderedNames()).toEqual(["B", "A", "C"]);
      fireEvent.click(screen.getAllByTestId("flag-button")[1]!); // A
      expect(screen.getByText(/2 values ready to commit · 1 skipped/i)).toBeInTheDocument();
      // First ⌘Z steps the SKIP back (LIFO).
      fireEvent.keyDown(window, { key: "z", metaKey: true });
      expect(screen.getByText(/3 values ready to commit/i)).toBeInTheDocument();
      expect(renderedNames()).toEqual(["B", "A", "C"]); // order still reordered
      // Second ⌘Z steps the REORDER back.
      fireEvent.keyDown(window, { key: "z", metaKey: true });
      expect(renderedNames()).toEqual(["A", "B", "C"]);
    });
  });
});
