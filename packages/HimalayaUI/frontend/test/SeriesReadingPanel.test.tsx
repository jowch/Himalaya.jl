/**
 * SeriesReadingPanel (Plan E, Task E-4 readout) — renders the derived
 * "phases present" reading: per-phase span + lattice trend rows, plus
 * coexistence / form-factor lines. Phase names are phase-coloured so
 * coexistence rows self-decode.
 */
import { describe, it, expect } from "vitest";
import { render } from "@testing-library/react";
import { SeriesReadingPanel } from "../src/components/SeriesReadingPanel";
import type { SeriesReading } from "../src/lib/series/seriesReading";

const reading: SeriesReading = {
  phases: [
    { phase: "Pn3m", spanLabel: "1:0 → 1:0.5", latticeTrend: "a 205 → 195 Å" },
    { phase: "Lamellar", spanLabel: "1:0.5 → 1:1", latticeTrend: "d 60 → 55 Å" },
  ],
  coexistenceAt: ["1:0.5"],
  formFactorOnlyAt: ["1:1.5"],
};

describe("SeriesReadingPanel", () => {
  it("renders one row per phase with its span + lattice trend", () => {
    const { getAllByTestId, getByText } = render(<SeriesReadingPanel reading={reading} />);
    expect(getAllByTestId("reading-phase-row")).toHaveLength(2);
    expect(getByText("a 205 → 195 Å")).toBeTruthy();
    expect(getByText("d 60 → 55 Å")).toBeTruthy();
  });

  it("renders coexistence and form-factor lines", () => {
    const { getByTestId } = render(<SeriesReadingPanel reading={reading} />);
    expect(getByTestId("reading-coexistence").textContent).toContain("1:0.5");
    expect(getByTestId("reading-form-factor").textContent).toContain("1:1.5");
  });

  it("phase names carry the phase colour via inline style", () => {
    const { getAllByTestId } = render(<SeriesReadingPanel reading={reading} />);
    const names = getAllByTestId("reading-phase-name");
    // The Pn3m name has a non-empty colour style.
    expect(names[0]!.getAttribute("style")).toMatch(/color/);
  });

  it("renders nothing for an empty reading", () => {
    const empty: SeriesReading = { phases: [], coexistenceAt: [], formFactorOnlyAt: [] };
    const { container } = render(<SeriesReadingPanel reading={empty} />);
    expect(container.querySelector('[data-testid="reading-phase-row"]')).toBeNull();
  });
});
