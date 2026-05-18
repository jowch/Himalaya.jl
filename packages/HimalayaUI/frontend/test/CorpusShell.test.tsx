import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import { CorpusShell } from "../src/components/CorpusShell";
import { SamplesPage } from "../src/pages/SamplesPage";

/** SamplesPage and CorpusTopbar both issue queries — answer them all 200/[]. */
function mockEmptyApi(): void {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify([]), {
      status: 200,
      headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("CorpusShell", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("renders the topbar and the matched child route in its Outlet", () => {
    mockEmptyApi();
    render(
      <QueryClientProvider client={makeClient()}>
        <MemoryRouter initialEntries={["/samples"]}>
          <Routes>
            <Route element={<CorpusShell />}>
              <Route path="/samples" element={<SamplesPage />} />
            </Route>
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );
    expect(screen.getByTestId("corpus-shell")).toBeInTheDocument();
    expect(screen.getByTestId("corpus-topbar")).toBeInTheDocument();
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
  });
});

describe("SamplesPage", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("renders the contact-sheet body", () => {
    mockEmptyApi();
    render(
      <QueryClientProvider client={makeClient()}>
        <MemoryRouter initialEntries={["/samples"]}>
          <SamplesPage />
        </MemoryRouter>
      </QueryClientProvider>,
    );
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
  });
});
