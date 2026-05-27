import { describe, it, expect, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter, useLocation } from "react-router-dom";
import type { ReactNode } from "react";
import { StaleUrlPage } from "../src/components/StaleUrlPage";
import { useAppState } from "../src/state";
import type { StaleUrlContext } from "../src/state";
import { _resetEmitMode } from "../src/lib/url/emitMode";

// StaleUrlPage uses useNavigate (its unknown-path CTA navigates to /samples
// since the Index surface was retired in #181), so renders need a Router.
let lastLocation = "";
function LocationProbe(): null {
  lastLocation = useLocation().pathname;
  return null;
}
function renderStale(ctx: StaleUrlContext, initialPath = "/foo/bar"): ReturnType<typeof render> {
  const wrap = (children: ReactNode) => (
    <MemoryRouter initialEntries={[initialPath]}>
      {children}
      <LocationProbe />
    </MemoryRouter>
  );
  return render(wrap(<StaleUrlPage staleUrlContext={ctx} />));
}

describe("StaleUrlPage", () => {
  beforeEach(() => {
    _resetEmitMode();
    lastLocation = "";
    useAppState.setState({
      staleUrlContext: null,
      activeExperimentId: undefined,
      activeSampleId: undefined,
      navModalOpen: false,
      navModalStep: "experiment",
    });
  });

  it("not_found:experiment — renders header, data-missing, CTA", () => {
    const ctx: StaleUrlContext = {
      kind: "not_found", missing: "experiment", missing_value: "lipid-typo",
      experiment_resolved: undefined, sample_resolved: undefined,
    };
    renderStale(ctx);
    expect(screen.getByTestId("stale-url-page")).toHaveAttribute("data-missing", "experiment");
    expect(screen.getByRole("heading")).toHaveTextContent(/Experiment 'lipid-typo' not found\./);
    fireEvent.click(screen.getByTestId("stale-url-cta"));
    expect(useAppState.getState().navModalOpen).toBe(true);
    expect(useAppState.getState().navModalStep).toBe("experiment");
    expect(useAppState.getState().staleUrlContext).toBeNull();
  });

  it("not_found:sample — preselects experiment and opens modal at sample step", () => {
    const ctx: StaleUrlContext = {
      kind: "not_found", missing: "sample", missing_value: "JC001",
      experiment_resolved: { id: 17, name: "lipid-screen" },
      sample_resolved: undefined,
    };
    renderStale(ctx);
    expect(screen.getByTestId("stale-url-page")).toHaveAttribute("data-missing", "sample");
    expect(screen.getByRole("heading")).toHaveTextContent(/Sample 'JC001' not found in 'lipid-screen'\./);
    fireEvent.click(screen.getByTestId("stale-url-cta"));
    const s = useAppState.getState();
    expect(s.activeExperimentId).toBe(17);
    expect(s.navModalOpen).toBe(true);
    expect(s.navModalStep).toBe("sample");
  });

  it("not_found:exposure — preserves sample, openModal=false (snap back)", () => {
    const ctx: StaleUrlContext = {
      kind: "not_found", missing: "exposure", missing_value: "JC001-099",
      experiment_resolved: { id: 17, name: "lipid-screen" },
      sample_resolved: { id: 42, name: "JC001" },
    };
    renderStale(ctx);
    expect(screen.getByTestId("stale-url-page")).toHaveAttribute("data-missing", "exposure");
    fireEvent.click(screen.getByTestId("stale-url-cta"));
    const s = useAppState.getState();
    expect(s.activeExperimentId).toBe(17);
    expect(s.activeSampleId).toBe(42);
    expect(s.navModalOpen).toBe(false);
  });

  it("unknown_path — Page not found, CTA navigates to /samples (#181)", () => {
    const ctx: StaleUrlContext = { kind: "unknown_path", raw: "/foo/bar" };
    renderStale(ctx);
    expect(screen.getByTestId("stale-url-page")).toHaveAttribute("data-missing", "path");
    expect(screen.getByRole("heading")).toHaveTextContent(/Page not found\./);
    expect(screen.getByTestId("stale-url-cta")).toHaveTextContent("Go to Samples");
    fireEvent.click(screen.getByTestId("stale-url-cta"));
    expect(lastLocation).toBe("/samples");
  });

  it("/ keypress triggers CTA", () => {
    const ctx: StaleUrlContext = {
      kind: "not_found", missing: "experiment", missing_value: "x",
      experiment_resolved: undefined, sample_resolved: undefined,
    };
    renderStale(ctx);
    fireEvent.keyDown(window, { key: "/" });
    expect(useAppState.getState().navModalOpen).toBe(true);
  });

  it("/ keypress inside an INPUT does NOT trigger CTA", () => {
    const ctx: StaleUrlContext = {
      kind: "not_found", missing: "experiment", missing_value: "x",
      experiment_resolved: undefined, sample_resolved: undefined,
    };
    renderStale(ctx);
    // Mount an input and dispatch keydown FROM that input as target.
    const input = document.createElement("input");
    document.body.appendChild(input);
    input.focus();
    fireEvent.keyDown(input, { key: "/" });
    expect(useAppState.getState().navModalOpen).toBe(false);
    document.body.removeChild(input);
  });
});
