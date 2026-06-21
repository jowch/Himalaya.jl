import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import { AppRoutes } from "../src/print/shell/AppRoutes";
import * as api from "../src/api";

// I1.7 (#163): Inspect is retired. Old /inspect* deep-links must land on the
// experiments home (/experiments), not the deleted Inspect page.
// T3.2: /inspect/* now redirects to /experiments (was /samples — which itself
// now redirects to /experiments).
describe("/inspect* redirect", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    // Mock the listExperiments call so ExperimentsHomePage renders.
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
  });

  it("redirects /inspect/:experiment/:sample to /experiments (T3.2)", async () => {
    render(
      <QueryClientProvider client={makeClient()}>
        <MemoryRouter initialEntries={["/inspect/1/10"]}>
          <AppRoutes />
        </MemoryRouter>
      </QueryClientProvider>,
    );
    await waitFor(() =>
      expect(screen.getByText("Your experiments")).toBeInTheDocument(),
    );
  });

  it("redirects bare /inspect to /experiments (T3.2)", async () => {
    render(
      <QueryClientProvider client={makeClient()}>
        <MemoryRouter initialEntries={["/inspect"]}>
          <AppRoutes />
        </MemoryRouter>
      </QueryClientProvider>,
    );
    await waitFor(() =>
      expect(screen.getByText("Your experiments")).toBeInTheDocument(),
    );
  });
});
