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

  it("R3-N1 (#209): drops the base card-header clamp + bottom border in the focus variant", () => {
    // Round-2 finding R3-N1: when headerSlot is supplied (focus variant),
    // the kicker+serif+subline 3-row stack overflows the 56px clamp and the
    // mockup `.plate-head` has no bottom hairline. Variant class must opt
    // into mockup-parity padding (`.card-header--slotted`).
    renderWithProviders(
      <PlotCard headerSlot={<div data-testid="custom-header">Custom</div>} />,
    );
    const strip = screen.getByTestId("plot-stat-strip");
    expect(strip).toHaveAttribute("data-variant", "slotted");
    expect(strip.className).toContain("card-header--slotted");
  });

  it("keeps the base card-header (no slotted variant) on the Index page", () => {
    // Sibling guardrail: prop-less PlotCard MUST keep the base card-header to
    // preserve top-edge alignment with the PhasePanel header on the Index page.
    renderWithProviders(<PlotCard />);
    const strip = screen.getByTestId("plot-stat-strip");
    expect(strip).toHaveAttribute("data-variant", "default");
    expect(strip.className).not.toContain("card-header--slotted");
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
