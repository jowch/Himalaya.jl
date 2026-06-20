import { describe, it, expect, vi, afterEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { CorpusTopbar } from "../src/print/shell/CorpusTopbar";
import * as api from "../src/api";

function at(path: string) {
  vi.spyOn(api, "listExperiments").mockResolvedValue([]);
  vi.spyOn(api, "listCorpusSamples").mockResolvedValue([]);
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}><CorpusTopbar /></MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("CorpusTopbar Experiments stage (Phase E1)", () => {
  afterEach(() => vi.restoreAllMocks());

  it("includes an Experiments stage tab linking to /experiments", () => {
    at("/samples");
    expect(screen.getByTestId("stage-tab-experiments")).toHaveAttribute("href", "/experiments");
  });
});
