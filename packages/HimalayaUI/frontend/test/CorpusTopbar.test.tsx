import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter, useSearchParams } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import { CorpusTopbar } from "../src/components/CorpusTopbar";

const EXPERIMENTS = [
  { id: 1, name: "SSRL Apr 2026", path: "/e1", data_dir: "/d1",
    analysis_dir: "/a1", manifest_path: null, created_at: "2026-04-01",
    q_units: "A-1" },
  { id: 2, name: "APS Jul 2026", path: "/e2", data_dir: "/d2",
    analysis_dir: "/a2", manifest_path: null, created_at: "2026-07-01",
    q_units: "A-1" },
];

function mockExperiments(): void {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify(EXPERIMENTS), {
      status: 200,
      headers: { "Content-Type": "application/json" },
    }),
  );
}

/** Surfaces the live ?beamtime= value so tests can assert URL writes. */
function BeamtimeProbe(): JSX.Element {
  const [params] = useSearchParams();
  return <span data-testid="beamtime-probe">{params.get("beamtime") ?? ""}</span>;
}

function renderTopbar(initialPath = "/samples") {
  return render(
    <QueryClientProvider client={makeClient()}>
      <MemoryRouter initialEntries={[initialPath]}>
        <CorpusTopbar />
        <BeamtimeProbe />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("CorpusTopbar", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("shows the corpus wordmark", () => {
    mockExperiments();
    renderTopbar();
    const wordmark = screen.getByTestId("corpus-wordmark");
    expect(wordmark).toHaveTextContent("Himalaya");
    expect(wordmark).toHaveTextContent("SAXS");
  });

  it("links the wordmark home to the corpus from any route", () => {
    mockExperiments();
    renderTopbar("/series");
    const wordmark = screen.getByTestId("corpus-wordmark");
    expect(wordmark.tagName).toBe("A");
    expect(wordmark).toHaveAttribute("href", "/samples");
  });

  it("renders three stage-tabs; Samples is active and links to /samples", () => {
    mockExperiments();
    renderTopbar();
    const samples = screen.getByTestId("stage-tab-samples");
    expect(samples).toHaveAttribute("href", "/samples");
    expect(samples).toHaveAttribute("data-active", "true");
    // Index remains inert until Phase 4; Series is now a live link (#173 / I3.3)
    // but inactive on /samples.
    expect(screen.getByTestId("stage-tab-index")).toBeDisabled();
    const series = screen.getByTestId("stage-tab-series");
    expect(series).toHaveAttribute("href", "/series");
    expect(series).not.toHaveAttribute("data-active");
  });

  it("lists experiments in the Beamtime chip once they load", async () => {
    mockExperiments();
    renderTopbar();
    await waitFor(() =>
      expect(screen.getByRole("option", { name: "SSRL Apr 2026" })).toBeInTheDocument(),
    );
    expect(screen.getByRole("option", { name: "APS Jul 2026" })).toBeInTheDocument();
  });

  it("writes ?beamtime= to the URL when an experiment is picked", async () => {
    mockExperiments();
    renderTopbar();
    await screen.findByRole("option", { name: "APS Jul 2026" });
    fireEvent.change(screen.getByTestId("beamtime-chip"), {
      target: { value: "2" },
    });
    expect(screen.getByTestId("beamtime-probe")).toHaveTextContent("2");
  });

  it("clears ?beamtime= when 'all' is picked", async () => {
    mockExperiments();
    renderTopbar("/samples?beamtime=2");
    await screen.findByRole("option", { name: "APS Jul 2026" });
    expect(screen.getByTestId("beamtime-probe")).toHaveTextContent("2");
    fireEvent.change(screen.getByTestId("beamtime-chip"), {
      target: { value: "" },
    });
    expect(screen.getByTestId("beamtime-probe")).toHaveTextContent("");
  });

  it("shows a Contact sheet | Loupe segmented switch", () => {
    mockExperiments();
    renderTopbar("/samples");
    const seg = screen.getByTestId("view-seg");
    expect(seg).toHaveTextContent("Contact sheet");
    expect(seg).toHaveTextContent("Loupe");
  });

  it("marks Contact sheet active on /samples and disables Loupe", () => {
    mockExperiments();
    renderTopbar("/samples");
    expect(screen.getByTestId("view-seg-sheet")).toHaveAttribute(
      "data-active",
      "true",
    );
    expect(screen.getByTestId("view-seg-loupe")).toBeDisabled();
  });

  it("marks Loupe active on a loupe route and links sheet back to /samples", () => {
    mockExperiments();
    renderTopbar("/samples/loupe/2");
    expect(screen.getByTestId("view-seg-loupe")).toHaveAttribute(
      "data-active",
      "true",
    );
    expect(screen.getByTestId("view-seg-sheet")).toHaveAttribute(
      "href",
      "/samples",
    );
  });

  it("preserves ?beamtime= on the Contact sheet link", () => {
    mockExperiments();
    renderTopbar("/samples/loupe/2?beamtime=1");
    expect(screen.getByTestId("view-seg-sheet")).toHaveAttribute(
      "href",
      "/samples?beamtime=1",
    );
  });

  it("reflects the current ?beamtime= as the selected option", async () => {
    mockExperiments();
    renderTopbar("/samples?beamtime=1");
    await screen.findByRole("option", { name: "SSRL Apr 2026" });
    // getByTestId returns HTMLElement (not generic); cast to read `.value`.
    expect(
      (screen.getByTestId("beamtime-chip") as HTMLSelectElement).value,
    ).toBe("1");
  });
});
