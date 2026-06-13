import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import type { Exposure, CorpusSample } from "../../src/api";

const setStatusMutate = vi.fn();
const selectMutate = vi.fn();
const addTagMutate = vi.fn();
const removeTagMutate = vi.fn();
const editTagMutate = vi.fn();

const state = {
  samples: [] as CorpusSample[],
  exposures: [] as Exposure[],
  loading: false,
};

vi.mock("../../src/queries", () => ({
  useCorpusSamples: () => ({ data: state.samples, isLoading: state.loading }),
  useExposures: () => ({ data: state.exposures, isLoading: state.loading }),
  useSetExposureStatus: () => ({ mutate: setStatusMutate }),
  useSelectExposure: () => ({ mutate: selectMutate }),
  useAddCorpusSampleTag: () => ({ mutate: addTagMutate }),
  useRemoveCorpusSampleTag: () => ({ mutate: removeTagMutate }),
  useEditCorpusSampleTag: () => ({ mutate: editTagMutate }),
  useCorpusSampleTags: () => ({ data: [] }),
}));

// boneyard Skeleton: render children when not loading (avoid capture machinery in JSDOM).
vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children }: { children: React.ReactNode }) => <>{children}</>,
}));

// DetectorImage touches fetch / createImageBitmap / OffscreenCanvas (absent in
// JSDOM); mock it at both the barrel and the direct path — LoupePage's own
// behaviour does not depend on its render.
vi.mock("../../src/print/detector", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));
vi.mock("../../src/print/detector/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));

import { LoupePage } from "../../src/print/pages/LoupePage";
import { setToastImpl } from "../../src/lib/toast";
import { setAnnounceImpl } from "../../src/lib/announce";

function exp(over: Partial<Exposure>): Exposure {
  return {
    id: 1, sample_id: 1, filename: "JC000-001.dat", kind: "file",
    selected: false, status: "accepted", image_path: "/x.tif", image_version: "v9",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null, ...over,
  };
}

/** Surfaces the live location.search so tests can pin "the loupe never
 *  rewrites the URL while flipping frames" (SA-F2 read-once contract). */
function LocationProbe(): JSX.Element {
  const loc = useLocation();
  return (
    <div
      data-testid="loc-probe"
      data-search={loc.search}
      data-pathname={loc.pathname}
    />
  );
}

function renderAt(sampleId: number, search = "") {
  return render(
    <MemoryRouter initialEntries={[`/samples/loupe/${sampleId}${search}`]}>
      <Routes>
        <Route
          path="/samples/loupe/:sampleId"
          element={
            <>
              <LoupePage />
              <LocationProbe />
            </>
          }
        />
        <Route path="/samples" element={<div data-testid="sheet">sheet</div>} />
      </Routes>
    </MemoryRouter>,
  );
}

/** Render the loupe with a router-state sample order (LO-NEXT), as the contact
 *  sheet hands it over. */
function renderWithOrder(sampleId: number, order: number[], search = "") {
  return render(
    <MemoryRouter
      initialEntries={[
        { pathname: `/samples/loupe/${sampleId}`, search, state: { sampleOrder: order } },
      ]}
    >
      <Routes>
        <Route
          path="/samples/loupe/:sampleId"
          element={
            <>
              <LoupePage />
              <LocationProbe />
            </>
          }
        />
        <Route path="/samples" element={<div data-testid="sheet">sheet</div>} />
      </Routes>
    </MemoryRouter>,
  );
}

/** Per-call mutate callbacks ({ onSuccess, onError }) captured by the mocks.
 *  Tests drive confirmation/failure explicitly — the toast contract is
 *  "announce on CONFIRMATION", never optimistically at mutate() time. */
type MutateCbs = { onSuccess?: (r: unknown) => void; onError?: (e: unknown) => void };
function lastCbs(mock: ReturnType<typeof vi.fn>): MutateCbs {
  return (mock.mock.calls.at(-1)?.[1] ?? {}) as MutateCbs;
}

beforeEach(() => {
  vi.clearAllMocks();
  state.samples = [{
    id: 42, experiment_id: 1, name: "JC042", display_name: "JC042 — LL37",
    notes: null, q_units: "A-1",
    tags: [{ id: 100, key: "LL37", value: "", source: "user" }],
  }];
  state.exposures = [
    exp({ id: 1, selected: true }),
    exp({ id: 2, status: "rejected" }),
  ];
  state.loading = false;
});

describe("LoupePage", () => {
  it("renders the headline, frame, side panel and filmstrip", () => {
    renderAt(42);
    expect(screen.getByTestId("loupe-page")).toBeInTheDocument();
    expect(screen.getByRole("heading", { name: /JC042 — LL37/ })).toBeInTheDocument();
    expect(screen.getByTestId("big-frame")).toBeInTheDocument();
    expect(screen.getByTestId("loupe-side-panel")).toBeInTheDocument();
  });

  it("sample title renders as h1 — the document outline has a top-level anchor (LO-NOH1)", () => {
    renderAt(42);
    // WCAG 1.3.1 / 2.4.6: the loupe has no outer h1 from a shell, so PlateHeader
    // must be promoted to h1 so assistive technology can navigate to the page topic.
    expect(screen.getByRole("heading", { level: 1, name: /JC042 — LL37/ })).toBeInTheDocument();
  });

  it("opens on the representative exposure (id 1)", () => {
    renderAt(42);
    expect(screen.getByTestId("big-frame")).not.toHaveAttribute("data-rejected");
  });

  it("the header subtitle speaks 'frame', not 'exposure' (LO-TERM consistency)", () => {
    const { container } = renderAt(42);
    const sub = container.querySelector('[data-role="plate-subtitle"]');
    expect(sub?.textContent).toMatch(/frame \d+ of \d+/);
    expect(sub?.textContent).not.toMatch(/exposure/i);
  });

  it("?exposure opens the loupe AT that frame (SA-F2)", () => {
    // Exposure 2 is the rejected, NON-default frame — landing on it proves the
    // param won over the representative default.
    const { container } = renderAt(42, "?exposure=2");
    expect(screen.getByTestId("big-frame")).toHaveAttribute("data-rejected", "true");
    expect(container.querySelector('[data-role="frame-caption"]')).toHaveTextContent(
      "frame 2 of 2",
    );
  });

  it("an unknown ?exposure falls back to the default frame silently", () => {
    // 999 is no frame of this sample. A FOREIGN sample's exposure id takes the
    // identical path: validation is membership in THIS sample's exposure list,
    // so foreign and unknown ids are indistinguishable here — both fall back.
    const { container } = renderAt(42, "?exposure=999");
    expect(screen.getByTestId("big-frame")).not.toHaveAttribute("data-rejected");
    expect(container.querySelector('[data-role="frame-caption"]')).toHaveTextContent(
      "frame 1 of 2",
    );
    // Silent fallback: no error surface of any kind.
    expect(screen.queryByTestId("loupe-not-found")).toBeNull();
  });

  it("re-selects a valid frame when the active one vanishes — never falsely 'no exposures' (LO-STALEACTIVEID)", () => {
    // Open on the representative (id 1), flip to the rejected frame (id 2), then
    // simulate another user dropping frame 2 via SSE (the exposure list loses
    // id 2). The loupe must heal to a valid frame, never strand on the gone id
    // and show the genuinely-frameless empty state while a frame still exists.
    // A FRESH element each render: a referentially-identical root element would
    // let React bail out of reconciling, so the SSE list change wouldn't apply.
    const tree = () => (
      <MemoryRouter initialEntries={["/samples/loupe/42"]}>
        <Routes>
          <Route path="/samples/loupe/:sampleId" element={<LoupePage />} />
        </Routes>
      </MemoryRouter>
    );
    const { rerender } = render(tree());
    fireEvent.keyDown(window, { key: "ArrowRight" }); // → frame 2 (rejected)
    expect(screen.getByTestId("big-frame")).toHaveAttribute("data-rejected", "true");

    // SSE drop: frame 2 is gone from this sample's exposures.
    state.exposures = [exp({ id: 1, selected: true })];
    rerender(tree());

    // No false empty state, and it healed to the surviving default frame.
    expect(screen.queryByText(/no exposures/i)).toBeNull();
    expect(screen.getByTestId("big-frame")).toBeInTheDocument();
    expect(screen.getByTestId("big-frame")).not.toHaveAttribute("data-rejected");
  });

  it("a malformed ?exposure is ignored (default frame, no error)", () => {
    const { container } = renderAt(42, "?exposure=abc");
    expect(container.querySelector('[data-role="frame-caption"]')).toHaveTextContent(
      "frame 1 of 2",
    );
  });

  it("frame flipping reads ?exposure once and never rewrites the URL", () => {
    renderAt(42, "?exposure=2");
    fireEvent.keyDown(window, { key: "ArrowLeft" }); // flip back to frame 1
    expect(screen.getByTestId("big-frame")).not.toHaveAttribute("data-rejected");
    // The param stays as the permalink wrote it — flipping is URL-silent.
    expect(screen.getByTestId("loc-probe")).toHaveAttribute(
      "data-search",
      "?exposure=2",
    );
  });

  it("drop toggle mutates status to rejected for the active exposure", () => {
    renderAt(42);
    fireEvent.keyDown(window, { key: "x" });
    expect(setStatusMutate).toHaveBeenCalledWith({ exposureId: 1, status: "rejected" }, expect.anything());
  });

  it("K marks an unscreened active frame accepted (SA-SCREENED)", () => {
    state.exposures = [exp({ id: 1, selected: true, status: null })];
    renderAt(42);
    fireEvent.keyDown(window, { key: "k" });
    expect(setStatusMutate).toHaveBeenCalledWith({ exposureId: 1, status: "accepted" }, expect.anything());
  });

  it("K on an accepted frame restores it to unscreened (toggle)", () => {
    renderAt(42); // fixture frame 1 is status "accepted"
    fireEvent.keyDown(window, { key: "k" });
    expect(setStatusMutate).toHaveBeenCalledWith({ exposureId: 1, status: null }, expect.anything());
  });

  it("K on a rejected frame sets accepted directly: last verb wins, no trip through null", () => {
    renderAt(42);
    fireEvent.keyDown(window, { key: "ArrowRight" }); // frame 2 is rejected
    fireEvent.keyDown(window, { key: "k" });
    expect(setStatusMutate).toHaveBeenCalledWith({ exposureId: 2, status: "accepted" }, expect.anything());
  });

  it("X on an accepted frame sets rejected directly: last verb wins", () => {
    renderAt(42); // fixture frame 1 is status "accepted"
    fireEvent.keyDown(window, { key: "x" });
    expect(setStatusMutate).toHaveBeenCalledWith({ exposureId: 1, status: "rejected" }, expect.anything());
  });

  it("K keep announces a toast ON CONFIRMATION, not before", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      state.exposures = [exp({ id: 1, selected: true, status: null })];
      renderAt(42);
      fireEvent.keyDown(window, { key: "k" });
      // No premature success claim while the save is still in flight.
      expect(toast).not.toHaveBeenCalled();
      lastCbs(setStatusMutate).onSuccess?.({});
      expect(toast).toHaveBeenCalledWith("Frame kept", "success");
    } finally {
      setToastImpl(null);
    }
  });

  it("K restore (accepted → null) announces 'Frame restored' on confirmation", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      renderAt(42); // frame 1 accepted
      fireEvent.keyDown(window, { key: "k" });
      lastCbs(setStatusMutate).onSuccess?.({});
      expect(toast).toHaveBeenCalledWith("Frame restored", "success");
    } finally {
      setToastImpl(null);
    }
  });

  it("dropping the LAST kept frame names the consequence + reversal, not a bland 'Frame dropped' (LO-LASTFRAME-NOGUARD)", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      // Default fixture: frame 1 accepted (the ONLY non-rejected frame), frame 2 rejected.
      renderAt(42);
      fireEvent.keyDown(window, { key: "x" }); // drop frame 1 — the last kept frame
      lastCbs(setStatusMutate).onSuccess?.({});
      expect(toast).toHaveBeenCalledWith(
        "Last kept frame dropped. This sample now has no kept frame; press X to restore it.",
        "warning",
      );
    } finally {
      setToastImpl(null);
    }
  });

  it("dropping a frame that is NOT the last kept one keeps the plain 'Frame dropped' toast", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      // Two kept frames — dropping one still leaves another, so no consequence copy.
      state.exposures = [
        exp({ id: 1, selected: true, status: "accepted" }),
        exp({ id: 2, status: "accepted" }),
      ];
      renderAt(42);
      fireEvent.keyDown(window, { key: "x" }); // drop frame 1; frame 2 stays kept
      lastCbs(setStatusMutate).onSuccess?.({});
      expect(toast).toHaveBeenCalledWith("Frame dropped", "success");
    } finally {
      setToastImpl(null);
    }
  });

  it("the Kept pill shows on an accepted frame; the Dropped pill on a rejected one", () => {
    const { container } = renderAt(42); // frame 1 accepted, frame 2 rejected
    expect(container.querySelector('[data-role="kept-tag"]')).toBeInTheDocument();
    expect(container.querySelector('[data-role="dropped-tag"]')).not.toBeInTheDocument();
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(container.querySelector('[data-role="kept-tag"]')).not.toBeInTheDocument();
    expect(container.querySelector('[data-role="dropped-tag"]')).toBeInTheDocument();
  });

  it("the frame caption reads the tri-state verdict word honestly", () => {
    state.exposures = [exp({ id: 1, selected: true, status: null })];
    const { container } = renderAt(42);
    expect(container.querySelector('[data-role="frame-caption"]')).toHaveTextContent(
      "unscreened",
    );
  });

  it("the loupe key legend documents K", () => {
    renderAt(42);
    expect(screen.getByText("keep / restore")).toBeInTheDocument();
  });

  it("R sets the representative when the active frame is NOT it", () => {
    renderAt(42);
    // Flip off the representative (frame 1) onto frame 2 first.
    fireEvent.keyDown(window, { key: "ArrowRight" });
    fireEvent.keyDown(window, { key: "r" });
    expect(selectMutate).toHaveBeenCalledWith(2, expect.anything());
  });

  it("R on the current representative is a no-op: no mutation, SR announce, no success toast (LO-REPLIES)", () => {
    const announce = vi.fn();
    const toast = vi.fn();
    setAnnounceImpl(announce);
    setToastImpl(toast);
    try {
      renderAt(42);
      // The loupe opens on the representative (exposure 1).
      fireEvent.keyDown(window, { key: "r" });
      expect(selectMutate).not.toHaveBeenCalled();
      expect(announce.mock.calls[0]?.[0]).toBe("Already the representative frame");
      expect(toast).not.toHaveBeenCalled();
    } finally {
      setAnnounceImpl(null);
      setToastImpl(null);
    }
  });

  it("ArrowRight flips to the next exposure", () => {
    renderAt(42);
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(screen.getByTestId("big-frame")).toHaveAttribute("data-rejected", "true");
  });

  it("X drop announces the success toast on CONFIRMATION, never before the save", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      // Two kept frames so this exercises the PLAIN drop path (not the
      // last-kept-frame consequence copy, covered by its own test).
      state.exposures = [
        exp({ id: 1, selected: true, status: "accepted" }),
        exp({ id: 2, status: "accepted" }),
      ];
      renderAt(42);
      fireEvent.keyDown(window, { key: "x" });
      // The save is still pending — claiming "Frame dropped" here would lie.
      expect(toast).not.toHaveBeenCalled();
      lastCbs(setStatusMutate).onSuccess?.({});
      expect(toast).toHaveBeenCalledWith("Frame dropped", "success");
    } finally {
      setToastImpl(null);
    }
  });

  it("X drop terminal infrastructure failure surfaces an error toast (no silent rollback)", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      renderAt(42);
      fireEvent.keyDown(window, { key: "x" });
      // No status on the error → network/5xx class, retries exhausted.
      lastCbs(setStatusMutate).onError?.(new TypeError("fetch failed"));
      expect(toast).toHaveBeenCalledTimes(1);
      expect(toast).toHaveBeenCalledWith(
        "Couldn't save the change. It was rolled back.",
        "error",
      );
    } finally {
      setToastImpl(null);
    }
  });

  it("X drop 4xx failure does NOT double-toast: the queue layer owns the validation toast", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      renderAt(42);
      fireEvent.keyDown(window, { key: "x" });
      lastCbs(setStatusMutate).onError?.(
        Object.assign(new Error("422"), { status: 422 }),
      );
      expect(toast).not.toHaveBeenCalled();
    } finally {
      setToastImpl(null);
    }
  });

  it("R set-representative announces a toast on confirmation (when it actually sets)", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      renderAt(42);
      fireEvent.keyDown(window, { key: "ArrowRight" });
      fireEvent.keyDown(window, { key: "r" });
      expect(toast).not.toHaveBeenCalled();
      lastCbs(selectMutate).onSuccess?.({});
      expect(toast).toHaveBeenCalledWith("Set as the representative frame", "success");
    } finally {
      setToastImpl(null);
    }
  });

  it("R set-representative terminal failure surfaces an error toast", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      renderAt(42);
      fireEvent.keyDown(window, { key: "ArrowRight" });
      fireEvent.keyDown(window, { key: "r" });
      lastCbs(selectMutate).onError?.(new TypeError("fetch failed"));
      expect(toast).toHaveBeenCalledWith(
        "Couldn't save the change. It was rolled back.",
        "error",
      );
    } finally {
      setToastImpl(null);
    }
  });

  it("the side panel omits the dead integration/collected placeholder rows (controls-don't-lie)", () => {
    renderAt(42);
    const panel = screen.getByTestId("loupe-side-panel");
    // The stubbed placeholder rows are never shown.
    expect(within(panel).queryByText("integration")).toBeNull();
    expect(within(panel).queryByText("collected")).toBeNull();
    // LO-EXPSPARSE: for a normal (non-rejected) frame the only meta row is the
    // position, which the subtitle + caption already carry, so the whole "This
    // frame" block is suppressed rather than rendering a lone redundant row.
    expect(within(panel).queryByText("This frame")).toBeNull();
  });

  it("a dropped frame with a rejection_reason tag shows the reason; other frames do not", () => {
    state.exposures = [
      exp({ id: 1, selected: true }),
      exp({
        id: 2,
        status: "rejected",
        tags: [{ id: 9, key: "rejection_reason", value: "beam flare", source: "manual" }],
      }),
    ];
    renderAt(42, "?exposure=2");
    const panel = screen.getByTestId("loupe-side-panel");
    expect(within(panel).getByText("reason")).toBeInTheDocument();
    expect(within(panel).getByText("beam flare")).toBeInTheDocument();
    // Flip to the kept frame — the reason row belongs to the dropped frame only.
    fireEvent.keyDown(window, { key: "ArrowLeft" });
    expect(within(panel).queryByText("reason")).toBeNull();
    expect(within(panel).queryByText("beam flare")).toBeNull();
  });

  it("a dropped frame WITHOUT a rejection_reason tag shows no reason row", () => {
    renderAt(42, "?exposure=2"); // fixture frame 2 is rejected, tagless
    const panel = screen.getByTestId("loupe-side-panel");
    expect(within(panel).queryByText("reason")).toBeNull();
  });

  it("rep-dropped sample shows the warning regardless of the active frame (LO-REPDROP)", () => {
    state.exposures = [
      exp({ id: 1, selected: true, status: "rejected" }),
      exp({ id: 2 }),
    ];
    renderAt(42);
    // Opens on the (dropped) representative.
    expect(screen.getByTestId("rep-dropped-warning")).toBeInTheDocument();
    // Still warned after flipping to a different, kept frame.
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(screen.getByTestId("rep-dropped-warning")).toBeInTheDocument();
  });

  it("no rep-dropped warning when the representative is kept", () => {
    renderAt(42);
    expect(screen.queryByTestId("rep-dropped-warning")).not.toBeInTheDocument();
    // A dropped NON-representative frame does not trigger it either.
    fireEvent.keyDown(window, { key: "ArrowRight" }); // frame 2 is rejected in the fixture
    expect(screen.queryByTestId("rep-dropped-warning")).not.toBeInTheDocument();
  });

  it("ArrowRight flip announces SR-only (frame position), not a toast", () => {
    const announce = vi.fn();
    const toast = vi.fn();
    setAnnounceImpl(announce);
    setToastImpl(toast);
    try {
      renderAt(42);
      fireEvent.keyDown(window, { key: "ArrowRight" });
      expect(announce.mock.calls[0]?.[0]).toBe("Frame 2 of 2");
      expect(toast).not.toHaveBeenCalled();
    } finally {
      setAnnounceImpl(null);
      setToastImpl(null);
    }
  });

  it("thumbnail-click navigation announces SR-only too (consistent with keyboard flip)", () => {
    const announce = vi.fn();
    const toast = vi.fn();
    setAnnounceImpl(announce);
    setToastImpl(toast);
    try {
      renderAt(42);
      // The filmstrip thumbnails; clicking the 2nd selects exposure 2.
      const thumbs = screen.getAllByTestId("thumbnail");
      fireEvent.click(thumbs[1]!);
      expect(announce.mock.calls[0]?.[0]).toBe("Frame 2 of 2");
      // Navigation is quiet (SR-only), never a toast.
      expect(toast).not.toHaveBeenCalled();
    } finally {
      setAnnounceImpl(null);
      setToastImpl(null);
    }
  });

  it("filmstrip marks kept on accepted frames ONLY (LO-KEPTSTRIP tri-state)", () => {
    state.exposures = [
      exp({ id: 1, selected: true, status: "accepted" }),
      exp({ id: 2, status: "rejected" }),
      exp({ id: 3, status: null }), // unscreened — neither kept nor dropped
    ];
    renderAt(42);
    expect(
      screen.getByRole("button", { name: "Frame 1, representative, kept" }),
    ).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "Frame 2, dropped" })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "Frame 3" })).toBeInTheDocument();
    const thumbs = screen.getAllByTestId("thumbnail");
    // The unscreened thumb stays "normal" — the regression this pin exists for.
    expect(thumbs[2]).toHaveAttribute("data-state", "normal");
  });

  it("Escape navigates back to the sheet", () => {
    renderAt(42);
    fireEvent.keyDown(window, { key: "Escape" });
    expect(screen.getByTestId("sheet")).toBeInTheDocument();
  });

  it("Back button navigates to the sheet", () => {
    renderAt(42);
    fireEvent.click(screen.getByTestId("loupe-back"));
    expect(screen.getByTestId("sheet")).toBeInTheDocument();
  });

  it("removing a sample tag calls removeTag with the tag's exact id (LO-TAGDUP)", () => {
    // beforeEach seeds sample 42 with one tag: id:100, key:"LL37", value:"" (key-only)
    renderAt(42);
    const panel = screen.getByTestId("loupe-side-panel");
    // The remove button for a key-only tag has aria-label "Remove LL37"
    const removeBtn = within(panel).getByRole("button", { name: /remove ll37/i });
    fireEvent.click(removeBtn);
    expect(removeTagMutate).toHaveBeenCalledWith(100);
  });

  it("removing the second of two byte-identical dose pills deletes tag id 7, not id 3 (LO-TAGDUP)", async () => {
    state.samples = [{
      id: 42, experiment_id: 1, name: "JC042", display_name: "JC042 — LL37",
      notes: null, q_units: "A-1",
      tags: [
        { id: 3, key: "dose", value: "10", source: "manual" },
        { id: 7, key: "dose", value: "10", source: "scoping" },
      ],
    }];
    renderAt(42);
    const panel = screen.getByTestId("loupe-side-panel");
    const removes = within(panel).getAllByRole("button", { name: /remove dose 10/i });
    // There are two "Remove dose 10" buttons — one per duplicate pill
    expect(removes).toHaveLength(2);
    // Click the second one (id 7, the scoping twin)
    fireEvent.click(removes[1]!);
    expect(removeTagMutate).toHaveBeenCalledWith(7);
    expect(removeTagMutate).not.toHaveBeenCalledWith(3);
  });

  it("shows not-found when the sample is missing", () => {
    state.samples = [];
    renderAt(999);
    expect(screen.getByTestId("loupe-not-found")).toBeInTheDocument();
  });

  it("not-found renders an EmptyState whose action goes back to the sheet (LP shell unification)", () => {
    state.samples = [];
    renderAt(999);
    const block = screen.getByTestId("loupe-not-found");
    expect(within(block).getByTestId("empty-state")).toBeInTheDocument();
    expect(
      within(block).getByRole("heading", { name: "Sample not found" }),
    ).toBeInTheDocument();
    fireEvent.click(within(block).getByRole("button", { name: "Back to the sheet" }));
    expect(screen.getByTestId("sheet")).toBeInTheDocument();
  });

  it("not-found keeps the Escape path back to the sheet", () => {
    state.samples = [];
    renderAt(999);
    fireEvent.keyDown(window, { key: "Escape" });
    expect(screen.getByTestId("sheet")).toBeInTheDocument();
  });

  it("shows the no-exposures state", () => {
    state.exposures = [];
    renderAt(42);
    expect(screen.getByText(/no exposures/i)).toBeInTheDocument();
  });

  it("modifier chords pass through: ⌘R does not set representative, ⌘X does not drop", () => {
    renderAt(42);
    // Flip OFF the representative first so a leaked plain-r WOULD mutate —
    // keeps this chord test load-bearing under the R no-op guard.
    fireEvent.keyDown(window, { key: "ArrowRight" });
    fireEvent.keyDown(window, { key: "r", metaKey: true });
    expect(selectMutate).not.toHaveBeenCalled();
    fireEvent.keyDown(window, { key: "x", metaKey: true });
    expect(setStatusMutate).not.toHaveBeenCalled();
    fireEvent.keyDown(window, { key: "x", ctrlKey: true });
    fireEvent.keyDown(window, { key: "r", altKey: true });
    // ⌘K belongs to the command palette / browser, never the keep verb.
    fireEvent.keyDown(window, { key: "k", metaKey: true });
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    expect(setStatusMutate).not.toHaveBeenCalled();
    expect(selectMutate).not.toHaveBeenCalled();
  });

  it("keys from a contenteditable target do not mutate", () => {
    renderAt(42);
    const editor = document.createElement("div");
    editor.setAttribute("contenteditable", "true");
    const inner = document.createElement("span");
    editor.appendChild(inner);
    document.body.appendChild(editor);
    try {
      fireEvent.keyDown(inner, { key: "x" });
      fireEvent.keyDown(inner, { key: "r" });
      fireEvent.keyDown(inner, { key: "k" });
      expect(setStatusMutate).not.toHaveBeenCalled();
      expect(selectMutate).not.toHaveBeenCalled();
    } finally {
      editor.remove();
    }
  });

  it("an open modal dialog suppresses X/R and Escape (no double action behind the modal)", () => {
    renderAt(42);
    const dialog = document.createElement("div");
    dialog.setAttribute("role", "dialog");
    dialog.setAttribute("aria-modal", "true");
    document.body.appendChild(dialog);
    try {
      fireEvent.keyDown(window, { key: "x" });
      fireEvent.keyDown(window, { key: "r" });
      expect(setStatusMutate).not.toHaveBeenCalled();
      expect(selectMutate).not.toHaveBeenCalled();
      // Escape must not ALSO navigate back while the modal owns it.
      fireEvent.keyDown(window, { key: "Escape" });
      expect(screen.queryByTestId("sheet")).toBeNull();
      expect(screen.getByTestId("loupe-page")).toBeInTheDocument();
    } finally {
      dialog.remove();
    }
  });
});

describe("LoupePage · LO-NEXT sample navigation", () => {
  const THREE = [
    { id: 10, experiment_id: 1, name: "S10", display_name: "Sample 10", notes: null, q_units: "A-1", tags: [] },
    { id: 11, experiment_id: 1, name: "S11", display_name: "Sample 11", notes: null, q_units: "A-1", tags: [] },
    { id: 12, experiment_id: 1, name: "S12", display_name: "Sample 12", notes: null, q_units: "A-1", tags: [] },
  ];
  beforeEach(() => {
    state.samples = THREE;
    state.exposures = [exp({ id: 1, selected: true })];
  });

  it("shows the position + enabled prev/next for a >1 sample walk handed over by the sheet", () => {
    renderWithOrder(11, [10, 11, 12]);
    expect(screen.getByTestId("loupe-sample-nav")).toBeInTheDocument();
    expect(screen.getByTestId("loupe-sample-position")).toHaveTextContent("2 / 3");
    expect(screen.getByTestId("loupe-prev-sample")).not.toBeDisabled();
    expect(screen.getByTestId("loupe-next-sample")).not.toBeDisabled();
  });

  it("Next walks to the next sample in the handed-over order", () => {
    renderWithOrder(11, [10, 11, 12]);
    fireEvent.click(screen.getByTestId("loupe-next-sample"));
    expect(screen.getByTestId("loc-probe")).toHaveAttribute("data-pathname", "/samples/loupe/12");
  });

  it("Prev walks to the previous sample in the handed-over order", () => {
    renderWithOrder(11, [10, 11, 12]);
    fireEvent.click(screen.getByTestId("loupe-prev-sample"));
    expect(screen.getByTestId("loc-probe")).toHaveAttribute("data-pathname", "/samples/loupe/10");
  });

  it("disables Prev on the first sample", () => {
    renderWithOrder(10, [10, 11, 12]);
    expect(screen.getByTestId("loupe-prev-sample")).toBeDisabled();
    expect(screen.getByTestId("loupe-next-sample")).not.toBeDisabled();
  });

  it("disables Next on the last sample", () => {
    renderWithOrder(12, [10, 11, 12]);
    expect(screen.getByTestId("loupe-next-sample")).toBeDisabled();
    expect(screen.getByTestId("loupe-prev-sample")).not.toBeDisabled();
  });

  it("']' steps to the next sample, '[' steps to the previous (arrows still flip frames)", () => {
    renderWithOrder(11, [10, 11, 12]);
    fireEvent.keyDown(window, { key: "]" });
    expect(screen.getByTestId("loc-probe")).toHaveAttribute("data-pathname", "/samples/loupe/12");
  });

  it("'[' steps to the previous sample", () => {
    renderWithOrder(11, [10, 11, 12]);
    fireEvent.keyDown(window, { key: "[" });
    expect(screen.getByTestId("loc-probe")).toHaveAttribute("data-pathname", "/samples/loupe/10");
  });

  it("renders no sample-nav when the walk holds a single sample", () => {
    renderWithOrder(11, [11]);
    expect(screen.queryByTestId("loupe-sample-nav")).not.toBeInTheDocument();
  });

  it("falls back to the beamtime-scoped corpus order on a direct URL (no router state)", () => {
    renderAt(11);
    expect(screen.getByTestId("loupe-sample-position")).toHaveTextContent("2 / 3");
    fireEvent.click(screen.getByTestId("loupe-next-sample"));
    expect(screen.getByTestId("loc-probe")).toHaveAttribute("data-pathname", "/samples/loupe/12");
  });

  it("documents the [ ] sample-step keys in the loupe legend", () => {
    renderWithOrder(11, [10, 11, 12]);
    expect(screen.getByText("prev / next sample")).toBeInTheDocument();
  });
});
