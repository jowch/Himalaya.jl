import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
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
