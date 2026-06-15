import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import { CorpusTopbar } from "../src/print/shell/CorpusTopbar";

function renderAt(path: string) {
  render(
    <QueryClientProvider client={makeClient()}>
      <MemoryRouter initialEntries={[path]}>
        <CorpusTopbar />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("CorpusTopbar — series tab", () => {
  it("renders the Series tab as a link to /series", () => {
    renderAt("/samples");
    const tab = screen.getByTestId("stage-tab-series");
    expect(tab).toHaveAttribute("href", "/series");
  });

  it("marks the Series tab active on a /series route", () => {
    renderAt("/series");
    expect(screen.getByTestId("stage-tab-series")).toHaveAttribute("aria-current", "page");
  });

  it("does not mark the Series tab active when on /samples", () => {
    renderAt("/samples");
    expect(screen.getByTestId("stage-tab-series")).not.toHaveAttribute("aria-current");
    // and the Samples tab IS active there
    expect(screen.getByTestId("stage-tab-samples")).toHaveAttribute("aria-current", "page");
  });
});
