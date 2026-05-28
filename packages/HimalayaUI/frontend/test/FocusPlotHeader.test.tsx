import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { FocusPlotHeader } from "../src/components/FocusPlotHeader";

describe("FocusPlotHeader", () => {
  it("renders the Integration kicker, serif sample name, and full provenance subline", () => {
    render(
      <FocusPlotHeader
        sampleName="Lipid 1-2 + LL37 1:0.5"
        sampleCode="smp_09"
        beamtime="SSRL Apr 2026"
        exposureLabel="smp_09_e03"
      />,
    );

    expect(screen.getByTestId("focus-plot-kicker")).toHaveTextContent("Integration");
    expect(screen.getByTestId("focus-plot-title")).toHaveTextContent(
      "Lipid 1-2 + LL37 1:0.5",
    );

    const sub = screen.getByTestId("focus-plot-sub").textContent ?? "";
    expect(sub).toContain("smp_09");
    expect(sub).toContain("SSRL Apr 2026");
    expect(sub).toContain("representative exposure smp_09_e03");
    // Segments joined with a middot separator.
    expect(sub).toContain("·");
  });

  it("never renders the experiment-picker placeholder text", () => {
    render(
      <FocusPlotHeader
        sampleName="C1"
        sampleCode="C1"
        beamtime="Beamtime A"
        exposureLabel="C1_e01"
      />,
    );
    expect(screen.queryByText(/pick an experiment/i)).toBeNull();
    expect(screen.queryByText(/click to change/i)).toBeNull();
  });

  it("R3-N2 (#209): subline max-width is 80ch, not the legacy 60ch", () => {
    // Round-2 finding R3-N2: 60ch truncated the most diagnostic suffix
    // ("representative exposure smp_NN_eNN") mid-token. The fix raises the
    // truncate cap to 80ch so the load-bearing exposure filename survives.
    render(
      <FocusPlotHeader
        sampleName="Lipid 1-2 + LL37 1:0.5"
        sampleCode="smp_09"
        beamtime="SSRL Apr 2026"
        exposureLabel="smp_09_e03"
      />,
    );
    const sub = screen.getByTestId("focus-plot-sub");
    expect(sub.className).toContain("max-w-[80ch]");
    expect(sub.className).not.toContain("max-w-[60ch]");
  });

  it("degrades to just the title when provenance fields are absent", () => {
    render(
      <FocusPlotHeader
        sampleName="C1"
        sampleCode={null}
        beamtime={null}
        exposureLabel={null}
      />,
    );
    expect(screen.getByTestId("focus-plot-title")).toHaveTextContent("C1");
    // No segments -> empty subline (no stray middots or "representative exposure").
    const sub = screen.getByTestId("focus-plot-sub").textContent ?? "";
    expect(sub.trim()).toBe("");
  });
});
