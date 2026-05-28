import { describe, it, expect, beforeEach, vi } from "vitest";
import { screen } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { PlotCard } from "../src/components/PlotCard";
import { useAppState } from "../src/state";

// PlotCard mounts queries (experiment/samples/trace/...) and reads app state;
// stub fetch so those settle to empty, and seed only the active sample so the
// card renders its "No trace data" body (header is what we assert on).
beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    activeSampleId: 1, activeExposureId: undefined,
    activeExperimentId: undefined, username: "tester",
  });
  vi.stubGlobal("ResizeObserver", class {
    observe() {} unobserve() {} disconnect() {}
  });
  vi.stubGlobal("fetch", vi.fn(async () =>
    new Response(JSON.stringify([]), {
      status: 200, headers: { "content-type": "application/json" },
    }),
  ));
});

describe("PlotCard headerSlot (focus variant)", () => {
  it("renders the supplied headerSlot in place of the experiment picker", () => {
    renderWithProviders(
      <PlotCard headerSlot={<div data-testid="custom-header">Custom</div>} />,
    );
    expect(screen.getByTestId("custom-header")).toBeInTheDocument();
    // The legacy experiment-picker button must NOT render in the focus variant.
    expect(screen.queryByTestId("plot-title")).toBeNull();
    expect(screen.queryByText(/pick an experiment/i)).toBeNull();
  });

  it("renders the experiment-picker title when no headerSlot is supplied", () => {
    renderWithProviders(<PlotCard />);
    expect(screen.getByTestId("plot-title")).toBeInTheDocument();
  });

  it("keeps the q-controls / export cluster in the focus variant", () => {
    renderWithProviders(
      <PlotCard headerSlot={<div data-testid="custom-header">Custom</div>} />,
    );
    // The x-scale toggle lives in the right-side cluster, outside the
    // header branch — it must remain regardless of the variant.
    expect(screen.getByTestId("x-scale-log")).toBeInTheDocument();
  });
});
