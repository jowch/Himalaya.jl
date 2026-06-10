import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import type { Exposure, CorpusSample } from "../../src/api";

const setStatusMutate = vi.fn();
const selectMutate = vi.fn();
const addTagMutate = vi.fn();
const removeTagMutate = vi.fn();

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

function renderAt(sampleId: number) {
  return render(
    <MemoryRouter initialEntries={[`/samples/loupe/${sampleId}`]}>
      <Routes>
        <Route path="/samples/loupe/:sampleId" element={<LoupePage />} />
        <Route path="/samples" element={<div data-testid="sheet">sheet</div>} />
      </Routes>
    </MemoryRouter>,
  );
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

  it("opens on the representative exposure (id 1)", () => {
    renderAt(42);
    expect(screen.getByTestId("big-frame")).not.toHaveAttribute("data-rejected");
  });

  it("drop toggle mutates status to rejected for the active exposure", () => {
    renderAt(42);
    fireEvent.keyDown(window, { key: "x" });
    expect(setStatusMutate).toHaveBeenCalledWith({ exposureId: 1, status: "rejected" });
  });

  it("R sets the representative", () => {
    renderAt(42);
    fireEvent.keyDown(window, { key: "r" });
    expect(selectMutate).toHaveBeenCalledWith(1);
  });

  it("ArrowRight flips to the next exposure", () => {
    renderAt(42);
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(screen.getByTestId("big-frame")).toHaveAttribute("data-rejected", "true");
  });

  it("X drop announces a status toast (consequential → visible)", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      renderAt(42);
      fireEvent.keyDown(window, { key: "x" });
      expect(toast).toHaveBeenCalledWith("Frame dropped", "success");
    } finally {
      setToastImpl(null);
    }
  });

  it("R set-representative announces a toast", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      renderAt(42);
      fireEvent.keyDown(window, { key: "r" });
      expect(toast).toHaveBeenCalledWith("Set as the representative frame", "success");
    } finally {
      setToastImpl(null);
    }
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

  it("removing a sample tag resolves its id via key+value", () => {
    renderAt(42);
    const panel = screen.getByTestId("loupe-side-panel");
    const removeBtn = within(panel).getByRole("button", { name: "Remove" });
    fireEvent.click(removeBtn);
    expect(removeTagMutate).toHaveBeenCalledWith(100);
  });

  it("shows not-found when the sample is missing", () => {
    state.samples = [];
    renderAt(999);
    expect(screen.getByTestId("loupe-not-found")).toBeInTheDocument();
  });

  it("shows the no-exposures state", () => {
    state.exposures = [];
    renderAt(42);
    expect(screen.getByText(/no exposures/i)).toBeInTheDocument();
  });

  it("modifier chords pass through: ⌘R does not set representative, ⌘X does not drop", () => {
    renderAt(42);
    fireEvent.keyDown(window, { key: "r", metaKey: true });
    expect(selectMutate).not.toHaveBeenCalled();
    fireEvent.keyDown(window, { key: "x", metaKey: true });
    expect(setStatusMutate).not.toHaveBeenCalled();
    fireEvent.keyDown(window, { key: "x", ctrlKey: true });
    fireEvent.keyDown(window, { key: "r", altKey: true });
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
