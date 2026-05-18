import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { CorpusTopbar } from "../src/components/CorpusTopbar";

function renderTopbar() {
  // MemoryRouter — the Samples stage-tab is a react-router <Link>.
  return render(
    <MemoryRouter>
      <CorpusTopbar />
    </MemoryRouter>,
  );
}

describe("CorpusTopbar", () => {
  it("shows the corpus wordmark", () => {
    renderTopbar();
    const wordmark = screen.getByTestId("corpus-wordmark");
    expect(wordmark).toHaveTextContent("Himalaya");
    expect(wordmark).toHaveTextContent("SAXS");
  });

  it("renders three stage-tabs", () => {
    renderTopbar();
    expect(screen.getByTestId("stage-tab-samples")).toBeInTheDocument();
    expect(screen.getByTestId("stage-tab-index")).toBeInTheDocument();
    expect(screen.getByTestId("stage-tab-series")).toBeInTheDocument();
  });

  it("makes Samples the active tab and links it to /samples", () => {
    renderTopbar();
    const samples = screen.getByTestId("stage-tab-samples");
    expect(samples).toHaveAttribute("href", "/samples");
    expect(samples).toHaveAttribute("data-active", "true");
  });

  it("renders Index and Series as inert (disabled) tabs", () => {
    renderTopbar();
    expect(screen.getByTestId("stage-tab-index")).toBeDisabled();
    expect(screen.getByTestId("stage-tab-series")).toBeDisabled();
  });

  it("shows the Beamtime facet chip", () => {
    renderTopbar();
    expect(screen.getByTestId("beamtime-chip")).toHaveTextContent("Beamtime");
  });
});
