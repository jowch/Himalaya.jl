import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import type { CorpusSample, Exposure } from "../src/api";
import { LoupePage } from "../src/pages/LoupePage";

// LoupePage is the only module under test that imports from ../src/queries.
// Mutable holders let each test set query state before render.
const h = vi.hoisted(() => ({
  corpusQ: {} as { data?: CorpusSample[]; isLoading: boolean },
  exposuresQ: {} as { data?: Exposure[]; isLoading: boolean },
  experimentQ: {} as { data?: { id: number; name: string | null; path: string } },
  setStatusMutate: vi.fn(),
  setRepMutate: vi.fn(),
}));

vi.mock("../src/queries", () => ({
  useCorpusSamples: () => h.corpusQ,
  useExposures: () => h.exposuresQ,
  useExperiment: () => h.experimentQ,
  useSetExposureStatus: () => ({
    mutate: h.setStatusMutate, isPending: false, error: null, reset: () => {},
  }),
  useSelectExposure: () => ({
    mutate: h.setRepMutate, isPending: false, error: null, reset: () => {},
  }),
}));

// DetectorImage touches fetch / createImageBitmap / OffscreenCanvas (absent in
// JSDOM); mock it — LoupePage's behaviour does not depend on its render.
vi.mock("../src/components/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));

function sample(over: Partial<CorpusSample> = {}): CorpusSample {
  return {
    id: 7, experiment_id: 3, name: "JC042", display_name: "DOPE 80%",
    notes: null, tags: [], q_units: "A-1", ...over,
  };
}

function exposure(over: Partial<Exposure> = {}): Exposure {
  return {
    id: 100, sample_id: 7, filename: "JC042-001.dat", kind: "file",
    selected: false, status: null, image_path: null, image_version: "",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
    ...over,
  };
}

function renderAt(path: string) {
  const client = makeClient();
  return render(
    <QueryClientProvider client={client}>
      <MemoryRouter initialEntries={[path]}>
        <Routes>
          <Route path="/samples" element={<div data-testid="samples-marker" />} />
          <Route path="/samples/loupe/:sampleId" element={<LoupePage />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("LoupePage — identity", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    h.corpusQ = { data: [sample()], isLoading: false };
    h.exposuresQ = { data: [exposure()], isLoading: false };
    h.experimentQ = { data: { id: 3, name: "Beamtime March", path: "/x" } };
  });

  it("renders the sample identified by the :sampleId route param", () => {
    renderAt("/samples/loupe/7");
    expect(screen.getByTestId("loupe-page")).toBeInTheDocument();
    expect(screen.getByText("DOPE 80%")).toBeInTheDocument();
    expect(screen.getByText(/Beamtime March/)).toBeInTheDocument();
  });
});

describe("LoupePage — composition", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    h.corpusQ = { data: [sample()], isLoading: false };
    h.exposuresQ = {
      data: [
        exposure({ id: 100, status: "accepted" }),
        exposure({ id: 101, status: "rejected" }),
      ],
      isLoading: false,
    };
    h.experimentQ = { data: { id: 3, name: "Beamtime March", path: "/x" } };
  });

  it("renders the frame and sidebar for a known sample", () => {
    renderAt("/samples/loupe/7");
    expect(screen.getByTestId("loupe-frame")).toBeInTheDocument();
    expect(screen.getByTestId("loupe-sidebar")).toBeInTheDocument();
  });

  it("default-selects the first accepted exposure", () => {
    renderAt("/samples/loupe/7");
    // Exposure 100 is accepted, 101 rejected → 100 is the default.
    expect(screen.getByTestId("loupe-meta-frame")).toHaveTextContent("1 of 2");
    expect(screen.getByTestId("loupe-meta-status")).toHaveTextContent("accepted");
  });

  it("flips the active exposure when a strip thumbnail is clicked", () => {
    renderAt("/samples/loupe/7");
    fireEvent.click(screen.getByTestId("thumb-cell-101"));
    expect(screen.getByTestId("loupe-meta-frame")).toHaveTextContent("2 of 2");
    expect(screen.getByTestId("loupe-meta-status")).toHaveTextContent("rejected");
  });

  it("shows a not-found panel for a sample id absent from the corpus", () => {
    renderAt("/samples/loupe/999");
    expect(screen.getByTestId("loupe-not-found")).toBeInTheDocument();
    expect(screen.queryByTestId("loupe-frame")).not.toBeInTheDocument();
  });

  it("navigates back to /samples when the back button is clicked", () => {
    renderAt("/samples/loupe/7");
    fireEvent.click(screen.getByTestId("loupe-back"));
    expect(screen.getByTestId("samples-marker")).toBeInTheDocument();
  });
});

describe("LoupePage — interactions", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    h.corpusQ = { data: [sample()], isLoading: false };
    h.exposuresQ = {
      data: [
        exposure({ id: 100, status: "accepted" }),
        exposure({ id: 101, status: null }),
      ],
      isLoading: false,
    };
    h.experimentQ = { data: { id: 3, name: "Beamtime March", path: "/x" } };
  });

  it("drops a kept exposure via the verdict toggle", () => {
    renderAt("/samples/loupe/7"); // default active = 100 (accepted)
    fireEvent.click(screen.getByTestId("loupe-drop-toggle"));
    expect(h.setStatusMutate).toHaveBeenCalledWith({
      exposureId: 100, status: "rejected",
    });
  });

  it("restores a dropped exposure via the verdict toggle", () => {
    h.exposuresQ = {
      data: [exposure({ id: 100, status: "rejected" })],
      isLoading: false,
    };
    renderAt("/samples/loupe/7");
    fireEvent.click(screen.getByTestId("loupe-drop-toggle"));
    expect(h.setStatusMutate).toHaveBeenCalledWith({
      exposureId: 100, status: null,
    });
  });

  it("sets the representative via the rep button", () => {
    renderAt("/samples/loupe/7");
    fireEvent.click(screen.getByTestId("loupe-set-representative"));
    expect(h.setRepMutate).toHaveBeenCalledWith(100);
  });

  it("flips frames with the arrow keys", () => {
    renderAt("/samples/loupe/7"); // default active = 100
    fireEvent.keyDown(document.body, { key: "ArrowRight" });
    expect(screen.getByTestId("loupe-meta-frame")).toHaveTextContent("2 of 2");
    fireEvent.keyDown(document.body, { key: "ArrowLeft" });
    expect(screen.getByTestId("loupe-meta-frame")).toHaveTextContent("1 of 2");
  });

  it("drops the active exposure with the X key", () => {
    renderAt("/samples/loupe/7");
    fireEvent.keyDown(document.body, { key: "x" });
    expect(h.setStatusMutate).toHaveBeenCalledWith({
      exposureId: 100, status: "rejected",
    });
  });

  it("sets the representative with the R key", () => {
    renderAt("/samples/loupe/7");
    fireEvent.keyDown(document.body, { key: "r" });
    expect(h.setRepMutate).toHaveBeenCalledWith(100);
  });

  it("goes back to /samples with the Escape key", () => {
    renderAt("/samples/loupe/7");
    fireEvent.keyDown(document.body, { key: "Escape" });
    expect(screen.getByTestId("samples-marker")).toBeInTheDocument();
  });
});

describe("LoupePage — loading", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    h.corpusQ = { isLoading: true };
    h.exposuresQ = { isLoading: true };
    h.experimentQ = {};
  });

  it("shows the loupe skeleton while data is loading", () => {
    renderAt("/samples/loupe/7");
    expect(screen.getByTestId("loupe-skeleton")).toBeInTheDocument();
    // Body content must not render while loading.
    expect(screen.queryByTestId("loupe-frame")).not.toBeInTheDocument();
    expect(screen.queryByTestId("loupe-not-found")).not.toBeInTheDocument();
  });
});
