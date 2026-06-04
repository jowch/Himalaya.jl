import { render, screen, fireEvent, within } from "@testing-library/react";
import { TracePlate } from "../../src/print/components/TracePlate";
import type { TraceModel } from "../../src/print/plot/TracePlot";

const model: TraceModel = {
  trace: { q: [0.02, 0.05, 0.1, 0.2], I: [10, 40, 20, 5], sigma: [1, 1, 1, 1] },
  peaks: [{ id: 0, q: 0.05, intensity: 40, source: "auto" }],
  phase: "Pn3m",
};

const base = {
  title: "Lipid 1-2 + LL37",
  trace: model,
  scale: "log" as const,
  onScaleChange: () => {},
};

describe("TracePlate", () => {
  it("renders the title and the trace plot region", () => {
    render(<TracePlate {...base} kicker="Integration" subtitle="smp_09" />);
    expect(screen.getByTestId("trace-plate")).toBeInTheDocument();
    expect(screen.getByText("Lipid 1-2 + LL37")).toBeInTheDocument();
    expect(screen.getByText("Integration")).toBeInTheDocument();
    expect(screen.getByText("smp_09")).toBeInTheDocument();
  });
  it("calls onScaleChange when the scale toggle is used", () => {
    const onScaleChange = vi.fn();
    render(<TracePlate {...base} onScaleChange={onScaleChange} />);
    fireEvent.click(screen.getByText("linear q"));
    expect(onScaleChange).toHaveBeenCalledWith("lin");
  });
  it("shows Auto-fit and fires onAutoFit", () => {
    const onAutoFit = vi.fn();
    render(<TracePlate {...base} onAutoFit={onAutoFit} />);
    fireEvent.click(screen.getByText("Auto-fit"));
    expect(onAutoFit).toHaveBeenCalled();
  });
  it("reflects the armed '+ Peak' toggle and fires onToggleAddPeak", () => {
    const onToggleAddPeak = vi.fn();
    render(<TracePlate {...base} addPeakArmed onToggleAddPeak={onToggleAddPeak} />);
    const peak = screen.getByText("+ Peak");
    expect(peak).toHaveAttribute("aria-pressed", "true");
    fireEvent.click(peak);
    expect(onToggleAddPeak).toHaveBeenCalled();
  });
  it("does NOT add a peak on a plot click when '+ Peak' is not armed", () => {
    // The arm gate: TracePlot would otherwise add a peak on any empty-space
    // click. With addPeakArmed=false the onAddPeak handler is withheld.
    const onAddPeak = vi.fn();
    const { container } = render(
      <TracePlate {...base} interaction={{ onXDomain: () => {}, onAddPeak }} />,
    );
    const svg = container.querySelector('svg[data-testid="trace-plate-plot"]')!;
    fireEvent.click(svg, { clientX: 300, clientY: 150 });
    expect(onAddPeak).not.toHaveBeenCalled();
  });
  it("adds a peak on a plot click only once armed", () => {
    const onAddPeak = vi.fn();
    const { container } = render(
      <TracePlate {...base} addPeakArmed interaction={{ onXDomain: () => {}, onAddPeak }} />,
    );
    const svg = container.querySelector('svg[data-testid="trace-plate-plot"]')!;
    fireEvent.click(svg, { clientX: 300, clientY: 150 });
    expect(onAddPeak).toHaveBeenCalledTimes(1);
    expect(typeof onAddPeak.mock.calls[0]![0]).toBe("number");
  });
  it("does NOT remove/exclude a peak on click when '+ Peak' is not armed", () => {
    // The arm governs all peak editing: while disarmed, clicking a peak must
    // not fire onClickPeak (remove / alt-exclude). Clicking the peak glyph at
    // q=0.05 (model fixture) would otherwise hit the remove path.
    const onClickPeak = vi.fn();
    const { container } = render(
      <TracePlate {...base} interaction={{ onXDomain: () => {}, onClickPeak }} />,
    );
    const glyph = container.querySelector('[data-role="plot-peaks"] [role="button"]');
    if (glyph) fireEvent.click(glyph);
    // Also drive a raw plot click in case the glyph isn't focusable in JSDOM —
    // either way no removal should fire while disarmed.
    const svg = container.querySelector('svg[data-testid="trace-plate-plot"]')!;
    fireEvent.click(svg, { clientX: 0, clientY: 0 });
    expect(onClickPeak).not.toHaveBeenCalled();
  });
  it("forwards a placement-only className", () => {
    render(<TracePlate {...base} className="mt-6" />);
    expect(screen.getByTestId("trace-plate").className).toContain("mt-6");
  });
  it("accepts a controlled xDomain (scroll-zoom window) and renders the plot", () => {
    // The zoom round-trip lives in the consumer: wheel → interaction.onXDomain →
    // store → xDomain. TracePlate's job is to plumb xDomain through to TracePlot;
    // the window itself is honoured by TracePlot (covered in its own suite).
    render(
      <TracePlate {...base} xDomain={[0.04, 0.06]} interaction={{ onXDomain: () => {} }} />,
    );
    expect(screen.getByTestId("trace-plate-plot")).toBeInTheDocument();
  });

  // q-link forwarding: an incoming hoveredQ near a peak (q=0.05) reaches the
  // inner TracePlot, which lights the matching peak → the q-readout chip appears.
  it("forwards hoveredQ to the inner TracePlot (q-readout chip lights up)", () => {
    const { container } = render(
      <TracePlate {...base} hoveredQ={0.05} interaction={{ onXDomain: () => {} }} />,
    );
    const readout = container.querySelector('[data-role="q-readout"]');
    expect(readout).toBeTruthy();
    expect(readout!.querySelector("text")!.textContent).toBe("0.050");
  });

  // onHoverQ forwarding: hovering (via the deterministic focus path) a peak glyph
  // bubbles the peak's q out through TracePlate's onHoverQ.
  it("forwards onHoverQ from the inner TracePlot", () => {
    const onHoverQ = vi.fn();
    const { container } = render(
      <TracePlate {...base} onHoverQ={onHoverQ} interaction={{ onXDomain: () => {} }} />,
    );
    const peakG = container.querySelector('[data-role="plot-peaks"] > g')!;
    onHoverQ.mockClear();
    fireEvent.focus(peakG);
    expect(onHoverQ).toHaveBeenCalledWith(0.05);
  });

  it("renders a node passed as actions inside the plate", () => {
    render(
      <TracePlate
        {...base}
        actions={<button data-testid="x" type="button">Export</button>}
      />,
    );
    const plate = screen.getByTestId("trace-plate");
    expect(within(plate).getByTestId("x")).toBeInTheDocument();
  });

  // highlightPeakIds forwarding: a peak NOT in the set dims (data-dimmed=true).
  it("forwards highlightPeakIds (a non-member peak dims)", () => {
    const twoPeaks: TraceModel = {
      ...model,
      peaks: [
        { id: 0, q: 0.05, intensity: 40, source: "auto" },
        { id: 1, q: 0.1, intensity: 20, source: "auto" },
      ],
    };
    const { container } = render(
      <TracePlate {...base} trace={twoPeaks} highlightPeakIds={new Set([0])} />,
    );
    const allGs = container.querySelectorAll('[data-role="plot-peaks"] > g');
    const g1 = Array.from(allGs).find((g) => g.querySelector('[data-peak-id="1"]'));
    expect(g1).toBeTruthy();
    expect(g1!.getAttribute("data-dimmed")).toBe("true");
  });
});
