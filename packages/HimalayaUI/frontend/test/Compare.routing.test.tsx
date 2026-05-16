/**
 * Compare.tsx — unified shell routing tests (Compare UX C-11).
 *
 * Verifies that all compare routes render the unified Compare shell
 * (identified by data-testid="compare-page"), and that the wrapper's
 * `display:contents` class preserves the underlying grid structure.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useAppState } from "../src/state";
import { AppShell } from "../src/components/AppShell";

function makeQc() {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

function renderShell(initialPath: string) {
  const qc = makeQc();
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[initialPath]}>
        <AppShell />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("Compare unified routing — Compare UX C-11", () => {
  beforeEach(() => {
    useAppState.setState({
      activePage: "compare",
      activeExperimentId: undefined,
      activeDraft: null,
    });
  });

  it("mounts Compare.tsx for /compare/all/:id", async () => {
    const { findByTestId } = renderShell("/compare/all/5");
    const el = await findByTestId("compare-page");
    expect(el).toBeInTheDocument();
  });

  it("mounts Compare.tsx for /compare/all/new", async () => {
    const { findByTestId } = renderShell("/compare/all/new");
    const el = await findByTestId("compare-page");
    expect(el).toBeInTheDocument();
  });

  it("mounts Compare.tsx for /experiments/:eid/compare/:id", async () => {
    useAppState.setState({ activePage: "compare", activeExperimentId: 3 });
    const { findByTestId } = renderShell("/experiments/3/compare/7");
    const el = await findByTestId("compare-page");
    expect(el).toBeInTheDocument();
  });

  it("mounts Compare.tsx for /experiments/:eid/compare/new", async () => {
    useAppState.setState({ activePage: "compare", activeExperimentId: 3 });
    const { findByTestId } = renderShell("/experiments/3/compare/new");
    const el = await findByTestId("compare-page");
    expect(el).toBeInTheDocument();
  });

  it("Compare.tsx wrapper uses display:contents so it does not take a layout slot", async () => {
    const { findByTestId } = renderShell("/compare/all/5");
    const el = await findByTestId("compare-page");
    // The wrapper must be `display:contents` — Tailwind class `contents`.
    // This ensures it does not consume a grid track in WorkspaceGrid.
    expect(el).toHaveClass("contents");
  });

  it("Compare.tsx delegates to review body when no draft is active", async () => {
    useAppState.setState({ activeDraft: null });
    const { findByTestId } = renderShell("/compare/all/5");
    await findByTestId("compare-page");
    // Review path: WorkspaceGrid is present inside the page body
    await waitFor(() => {
      expect(screen.getByTestId("workspace-grid")).toBeInTheDocument();
    });
  });
});
