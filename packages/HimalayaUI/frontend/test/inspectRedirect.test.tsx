import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import { AppRoutes } from "../src/components/AppRoutes";

// I1.7 (#163): Inspect is retired. Old /inspect* deep-links must land on the
// corpus contact sheet (/samples), not the deleted Inspect page.
describe("/inspect* redirect", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    // All corpus/experiments fetches return empty so SamplesPage reaches its
    // empty state rather than erroring — the redirect target is what matters.
    vi.spyOn(global, "fetch").mockImplementation(() =>
      Promise.resolve(
        new Response("[]", {
          status: 200,
          headers: { "Content-Type": "application/json" },
        }),
      ),
    );
  });

  it("redirects /inspect/:experiment/:sample to /samples", async () => {
    render(
      <QueryClientProvider client={makeClient()}>
        <MemoryRouter initialEntries={["/inspect/1/10"]}>
          <AppRoutes />
        </MemoryRouter>
      </QueryClientProvider>,
    );
    await waitFor(() =>
      expect(screen.getByTestId("samples-page")).toBeInTheDocument(),
    );
  });

  it("redirects bare /inspect to /samples", async () => {
    render(
      <QueryClientProvider client={makeClient()}>
        <MemoryRouter initialEntries={["/inspect"]}>
          <AppRoutes />
        </MemoryRouter>
      </QueryClientProvider>,
    );
    await waitFor(() =>
      expect(screen.getByTestId("samples-page")).toBeInTheDocument(),
    );
  });
});
