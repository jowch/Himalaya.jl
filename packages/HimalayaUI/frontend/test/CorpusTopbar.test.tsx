import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter, useSearchParams } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import { CorpusTopbar } from "../src/print/shell/CorpusTopbar";

const EXPERIMENTS = [
  { id: 1, name: "SSRL Apr 2026", path: "/e1", data_dir: "/d1",
    analysis_dir: "/a1", manifest_path: null, created_at: "2026-04-01",
    q_units: "A-1" },
  { id: 2, name: "APS Jul 2026", path: "/e2", data_dir: "/d2",
    analysis_dir: "/a2", manifest_path: null, created_at: "2026-07-01",
    q_units: "A-1" },
];

function mockExperiments(): void {
  // A fresh Response per call: the topbar runs TWO queries (experiments +
  // corpus samples, both needed by the SA-F5 unknown-filter verdict), and a
  // Response body is single-use — a shared mockResolvedValue would leave the
  // second query permanently unresolved ("body already used").
  vi.spyOn(global, "fetch").mockImplementation(() =>
    Promise.resolve(
      new Response(JSON.stringify(EXPERIMENTS), {
        status: 200,
        headers: { "Content-Type": "application/json" },
      }),
    ),
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

  it("renders two stage-tabs (Samples, Series); Samples is active and links to /samples", () => {
    mockExperiments();
    renderTopbar();
    const samples = screen.getByTestId("stage-tab-samples");
    expect(samples).toHaveAttribute("href", "/samples");
    expect(samples).toHaveAttribute("data-active", "true");
    // The focus workspace is no longer a top-level stage tab — it is reached by
    // opening a sample. Only Samples + Series remain.
    expect(screen.queryByTestId("stage-tab-index")).toBeNull();
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

  // The beamtime filter is honored only on the samples surface. It used to
  // render on every route — changeable but inert on /series and /sample/:id.
  // Hide it where it does nothing rather than lie.
  it("shows the beamtime chip on the samples surface", () => {
    mockExperiments();
    renderTopbar("/samples");
    expect(screen.getByTestId("beamtime-chip")).toBeInTheDocument();
  });

  it("hides the beamtime chip on /series (where it is inert)", () => {
    mockExperiments();
    renderTopbar("/series");
    expect(screen.queryByTestId("beamtime-chip")).toBeNull();
  });

  it("hides the beamtime chip on the focus route (where it is inert)", () => {
    mockExperiments();
    renderTopbar("/sample/7");
    expect(screen.queryByTestId("beamtime-chip")).toBeNull();
  });

  it("an unknown ?beamtime shows the honest disabled option, never 'all experiments' (SA-F5)", async () => {
    mockExperiments();
    renderTopbar("/samples?beamtime=99");
    // Wait for BOTH lists (experiments + corpus) to settle the unknown verdict.
    const unknown = await screen.findByRole("option", { name: "Unknown beamtime" });
    expect(unknown).toBeDisabled();
    const select = screen.getByTestId("beamtime-chip") as HTMLSelectElement;
    expect(select.value).toBe("99");
    // The displayed selection is the honest label — the same copy string the
    // contact sheet's EmptyState shows (shared UNKNOWN_BEAMTIME_LABEL).
    expect(select.selectedOptions[0]?.textContent).toBe("Unknown beamtime");
  });

  it("picking 'all experiments' is the way out of an unknown ?beamtime", async () => {
    mockExperiments();
    renderTopbar("/samples?beamtime=99");
    await screen.findByRole("option", { name: "Unknown beamtime" });
    fireEvent.change(screen.getByTestId("beamtime-chip"), { target: { value: "" } });
    expect(screen.getByTestId("beamtime-probe")).toHaveTextContent("");
    await waitFor(() =>
      expect(screen.queryByRole("option", { name: "Unknown beamtime" })).toBeNull(),
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
