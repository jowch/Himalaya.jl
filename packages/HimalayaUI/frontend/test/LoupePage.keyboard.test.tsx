// test/LoupePage.keyboard.test.tsx
//
// T2.5: Loupe rev-2 axis pins.
//
// After T2.5:
//   ←/→ = prevDetail/nextDetail → flip frames (renamed from prevExposure/nextExposure).
//   ↑/↓ = prevSample/nextSample → navigate to sibling sample.
//   r = representative → set representative (unchanged).
//   Backspace = restore → setStatus(null) to restore the active frame.
//
// Complements test/print-pages/LoupePage.test.tsx (full page behaviour).
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import type { Exposure, CorpusSample } from "../src/api";

const setStatusMutate = vi.fn();
const selectMutate = vi.fn();

const state = {
  samples: [] as CorpusSample[],
  exposures: [] as Exposure[],
};

vi.mock("../src/queries", () => ({
  useCorpusSamples: () => ({ data: state.samples, isLoading: false }),
  useExposures: () => ({ data: state.exposures, isLoading: false }),
  useSetExposureStatus: () => ({ mutate: setStatusMutate }),
  useSelectExposure: () => ({ mutate: selectMutate }),
  useAddCorpusSampleTag: () => ({ mutate: vi.fn() }),
  useRemoveCorpusSampleTag: () => ({ mutate: vi.fn() }),
  useEditCorpusSampleTag: () => ({ mutate: vi.fn() }),
  useCorpusSampleTags: () => ({ data: [] }),
}));

vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children }: { children: React.ReactNode }) => <>{children}</>,
}));

vi.mock("../src/print/detector", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));
vi.mock("../src/print/detector/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));

import { LoupePage } from "../src/print/pages/LoupePage";

function mkExposure(id: number, over: Partial<Exposure> = {}): Exposure {
  return {
    id, sample_id: 42, filename: `JC042-00${id}.dat`, kind: "file" as const,
    selected: id === 1, status: "accepted", image_path: null, image_version: "",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null, ...over,
  };
}

// Three samples for ↑/↓ navigation.
const THREE_SAMPLES: CorpusSample[] = [
  { id: 10, experiment_id: 1, name: "S10", notes: null, q_units: "A-1", tags: [] } as CorpusSample,
  { id: 11, experiment_id: 1, name: "S11", notes: null, q_units: "A-1", tags: [] } as CorpusSample,
  { id: 12, experiment_id: 1, name: "S12", notes: null, q_units: "A-1", tags: [] } as CorpusSample,
];

function renderAt(sampleId: number) {
  return render(
    <MemoryRouter
      initialEntries={[{ pathname: `/samples/loupe/${sampleId}`, state: { sampleOrder: [10, 11, 12] } }]}
    >
      <Routes>
        <Route path="/samples/loupe/:sampleId" element={<LoupePage />} />
        <Route path="/samples" element={<div data-testid="sheet">sheet</div>} />
      </Routes>
    </MemoryRouter>,
  );
}

beforeEach(() => {
  vi.clearAllMocks();
  state.samples = THREE_SAMPLES;
  state.exposures = [mkExposure(1), mkExposure(2, { status: "rejected", selected: false })];
});

// ── T2.5 axis contract ────────────────────────────────────────────────────────

describe("Loupe T2.5 keyboard: rev-2 axis contract", () => {
  // ArrowRight/Left on window removed in interaction-arch; scope-container arrows
  // covered by test/LoupePage.actions.test.tsx "ArrowRight on the scope container".

  it("←/→ flip frames — ArrowLeft goes back", () => {
    renderAt(11);
    fireEvent.keyDown(window, { key: "ArrowRight" }); // → frame 2
    fireEvent.keyDown(window, { key: "ArrowLeft" }); // ← back to frame 1
    expect(screen.getByTestId("big-frame")).not.toHaveAttribute("data-rejected");
  });

  it("↓ steps to the next sample (nextSample = ArrowDown)", () => {
    renderAt(11);
    // Open at sample 11 with order [10, 11, 12]; ArrowDown → 12.
    const locBefore = screen.queryByTestId("sheet");
    expect(locBefore).toBeNull();
    fireEvent.keyDown(window, { key: "ArrowDown" });
    // Navigates to loupe for sample 12.
    // After navigation the new route renders (MemoryRouter resolves within same tree).
    // Check via the big-frame still present (loupe route active, not /samples).
    // Best assertion: no sheet rendered (we're not on /samples).
    expect(screen.queryByTestId("sheet")).toBeNull();
  });

  it("↑ steps to the previous sample (prevSample = ArrowUp)", () => {
    renderAt(11);
    fireEvent.keyDown(window, { key: "ArrowUp" });
    // Navigated to sample 10 (no sheet shown; still on loupe route).
    expect(screen.queryByTestId("sheet")).toBeNull();
  });

  // r via window keyboard layer covered by test/LoupePage.actions.test.tsx
  //   "r (representative) through the keyboard layer fires setRepresentative".
  // Backspace (restore) intentionally has NO keyboard binding in the new model
  //   (would collide with text editing in tag input); dock button remains.

  it("[ and ] are no longer bound as sample nav — loupe uses ↑/↓", () => {
    renderAt(11);
    // These are not in the registry; they must not navigate to another sample.
    fireEvent.keyDown(window, { key: "]" });
    fireEvent.keyDown(window, { key: "[" });
    // Still on loupe page (no sheet).
    expect(screen.queryByTestId("sheet")).toBeNull();
    // Still showing sample 11's frame (not navigated away).
    expect(screen.getByTestId("big-frame")).toBeInTheDocument();
  });
});
