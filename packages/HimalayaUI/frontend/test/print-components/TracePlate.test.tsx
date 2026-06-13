import { render, screen, fireEvent, within } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
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

  it("names the trace figure svg so it is not a nameless img in the a11y tree (WCAG 1.1.1)", () => {
    const { container } = render(<TracePlate {...base} />);
    const svg = container.querySelector('svg[data-testid="trace-plate-plot"]')!;
    expect(svg.getAttribute("role")).toBe("img");
    expect(svg.getAttribute("aria-label")).toBe("Integration trace: intensity vs q");
  });

  // ── H10: disarmed discoverability cue ────────────────────────────────────────

  it("shows the disarmed discoverability cue when the edit toggle is wired and not armed", () => {
    render(<TracePlate {...base} onToggleAddPeak={() => {}} />);
    expect(screen.getByTestId("peak-edit-cue").textContent).toBe(
      "Arm + Peak to edit peaks.",
    );
  });

  it("hides the disarmed cue once armed (never both the cue and the armed legend)", () => {
    render(
      <TracePlate
        {...base}
        addPeakArmed
        onToggleAddPeak={() => {}}
        interaction={{ onXDomain: () => {}, onAddPeak: () => {}, onClickPeak: () => {} }}
      />,
    );
    expect(screen.queryByTestId("peak-edit-cue")).toBeNull();
    expect(screen.getByTestId("peak-edit-hint")).toBeInTheDocument();
  });

  it("omits the disarmed cue when no edit affordance exists (no onToggleAddPeak)", () => {
    render(<TracePlate {...base} />);
    expect(screen.queryByTestId("peak-edit-cue")).toBeNull();
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

  // ── FO-EDIT: discoverability hint + minimal keyboard path ────────────────────

  const HINT_RX = /Click the trace to add a peak/;

  it("hides the edit hint (and the add-at-q field) when '+ Peak' is not armed", () => {
    render(
      <TracePlate
        {...base}
        interaction={{ onXDomain: () => {}, onAddPeak: () => {}, onClickPeak: () => {} }}
      />,
    );
    expect(screen.queryByText(HINT_RX)).toBeNull();
    expect(screen.queryByTestId("add-peak-at-q")).toBeNull();
  });

  it("shows the three-verb edit hint while armed", () => {
    render(
      <TracePlate
        {...base}
        addPeakArmed
        interaction={{ onXDomain: () => {}, onAddPeak: () => {}, onClickPeak: () => {} }}
      />,
    );
    expect(
      screen.getByText(
        "Click the trace to add a peak. Click a peak to remove it. Alt-click excludes it from indexing.",
      ),
    ).toBeInTheDocument();
  });

  it("unarmed peak marks are read-only: no role, no tabindex", () => {
    const { container } = render(
      <TracePlate
        {...base}
        interaction={{ onXDomain: () => {}, onAddPeak: () => {}, onClickPeak: () => {} }}
      />,
    );
    expect(
      container.querySelector('[data-role="plot-peaks"] [role="button"]'),
    ).toBeNull();
    expect(
      container.querySelector('[data-role="plot-peaks"] [tabindex]'),
    ).toBeNull();
  });

  it("armed peak marks are focusable buttons named by their q", () => {
    const { container } = render(
      <TracePlate
        {...base}
        addPeakArmed
        interaction={{ onXDomain: () => {}, onAddPeak: () => {}, onClickPeak: () => {} }}
      />,
    );
    const mark = container.querySelector('[data-role="plot-peaks"] [role="button"]');
    expect(mark).toBeTruthy();
    expect(mark!.getAttribute("tabindex")).toBe("0");
    expect(mark!.getAttribute("aria-label")).toBe("Auto peak at q = 0.0500");
  });

  it("armed: Enter on a peak mark removes it; Alt+Enter excludes it", () => {
    const onClickPeak = vi.fn();
    const { container } = render(
      <TracePlate
        {...base}
        addPeakArmed
        interaction={{ onXDomain: () => {}, onClickPeak }}
      />,
    );
    const mark = container.querySelector('[data-role="plot-peaks"] [role="button"]')!;
    fireEvent.keyDown(mark, { key: "Enter" });
    expect(onClickPeak).toHaveBeenLastCalledWith(0, false);
    fireEvent.keyDown(mark, { key: "Enter", altKey: true });
    expect(onClickPeak).toHaveBeenLastCalledWith(0, true);
  });

  it("armed add-at-q: typing a q inside the trace domain and submitting fires onAddPeak(q)", () => {
    const onAddPeak = vi.fn();
    render(
      <TracePlate
        {...base}
        addPeakArmed
        interaction={{ onXDomain: () => {}, onAddPeak }}
      />,
    );
    const input = screen.getByLabelText("q value for new peak");
    fireEvent.change(input, { target: { value: "0.07" } });
    const add = screen.getByRole("button", { name: "Add peak at q" });
    expect(add).toBeEnabled();
    fireEvent.click(add);
    expect(onAddPeak).toHaveBeenCalledTimes(1);
    expect(onAddPeak).toHaveBeenCalledWith(0.07);
  });

  it("armed add-at-q: a q outside the trace domain disables Add and never fires", () => {
    // model trace q spans [0.02, 0.2]
    const onAddPeak = vi.fn();
    render(
      <TracePlate
        {...base}
        addPeakArmed
        interaction={{ onXDomain: () => {}, onAddPeak }}
      />,
    );
    const input = screen.getByLabelText("q value for new peak");
    const add = screen.getByRole("button", { name: "Add peak at q" });
    // Empty → disabled.
    expect(add).toBeDisabled();
    fireEvent.change(input, { target: { value: "0.5" } });
    expect(add).toBeDisabled();
    fireEvent.click(add);
    expect(onAddPeak).not.toHaveBeenCalled();
    fireEvent.change(input, { target: { value: "0.001" } });
    expect(add).toBeDisabled();
  });

  it("armed add-at-q while ZOOMED validates against the visible window, not the full extent (FO-ZOOMEDIT)", () => {
    // Trace q spans [0.02, 0.2]; zoom the visible window to [0.04, 0.06]. A q of
    // 0.1 is in the data but OFF-SCREEN — adding it would seed a peak invisible
    // until zoom-out. The field must only accept q within what is currently
    // shown.
    const onAddPeak = vi.fn();
    render(
      <TracePlate
        {...base}
        addPeakArmed
        xDomain={[0.04, 0.06]}
        interaction={{ onXDomain: () => {}, onAddPeak }}
      />,
    );
    const input = screen.getByLabelText("q value for new peak");
    const add = screen.getByRole("button", { name: "Add peak at q" });
    // In the data extent but outside the visible window → rejected.
    fireEvent.change(input, { target: { value: "0.1" } });
    expect(add).toBeDisabled();
    fireEvent.click(add);
    expect(onAddPeak).not.toHaveBeenCalled();
    // Inside the visible window → accepted.
    fireEvent.change(input, { target: { value: "0.05" } });
    expect(add).toBeEnabled();
    fireEvent.click(add);
    expect(onAddPeak).toHaveBeenCalledWith(0.05);
  });

  it("onClickPeak only: hint drops the add sentence and the add-at-q field", () => {
    // The hint promises only the wired verbs: with no onAddPeak, "Click the
    // trace to add a peak" would be false, and no field may promise an add
    // path the plate cannot deliver.
    render(
      <TracePlate
        {...base}
        addPeakArmed
        interaction={{ onXDomain: () => {}, onClickPeak: () => {} }}
      />,
    );
    const hint = screen.getByTestId("peak-edit-hint");
    expect(hint.textContent).toBe(
      "Click a peak to remove it. Alt-click excludes it from indexing.",
    );
    expect(screen.queryByText(HINT_RX)).toBeNull();
    expect(screen.queryByTestId("add-peak-at-q")).toBeNull();
  });

  it("onAddPeak only: hint drops the remove/exclude sentences", () => {
    // Worse than false: without onClickPeak wired, TracePlot's click handler
    // falls through to onAddPeak on a peak hit, so "Click a peak to remove it"
    // would actually duplicate-add. The hint must not promise it.
    render(
      <TracePlate
        {...base}
        addPeakArmed
        interaction={{ onXDomain: () => {}, onAddPeak: () => {} }}
      />,
    );
    const hint = screen.getByTestId("peak-edit-hint");
    expect(within(hint).getByText(HINT_RX)).toBeInTheDocument();
    expect(hint.textContent).not.toMatch(/remove|excludes/);
    expect(screen.getByTestId("add-peak-at-q")).toBeInTheDocument();
  });

  it("clears a typed-but-unsubmitted q when the plate disarms", () => {
    const interaction = { onXDomain: () => {}, onAddPeak: () => {} };
    const { rerender } = render(
      <TracePlate {...base} addPeakArmed interaction={interaction} />,
    );
    const input = screen.getByLabelText("q value for new peak");
    fireEvent.change(input, { target: { value: "0.07" } });
    expect(input).toHaveValue(0.07);
    // Disarm, then re-arm: the stale draft must not reappear.
    rerender(<TracePlate {...base} interaction={interaction} />);
    rerender(<TracePlate {...base} addPeakArmed interaction={interaction} />);
    expect(screen.getByLabelText("q value for new peak")).toHaveValue(null);
  });

  // ── F7: Escape disarms the armed + Peak mode ─────────────────────────────────
  // The armed mode had no keyboard exit. Escape disarms via onToggleAddPeak,
  // with a modal dialog winning (ModalShell owns Escape-to-close and stamps
  // preventDefault on the press it consumes), and the hint names the exit.

  it("Escape disarms: fires onToggleAddPeak while armed", () => {
    const onToggleAddPeak = vi.fn();
    render(
      <TracePlate
        {...base}
        addPeakArmed
        onToggleAddPeak={onToggleAddPeak}
        interaction={{ onXDomain: () => {}, onAddPeak: () => {} }}
      />,
    );
    fireEvent.keyDown(document.body, { key: "Escape" });
    expect(onToggleAddPeak).toHaveBeenCalledTimes(1);
  });

  it("Escape-disarm re-anchors focus to '+ Peak' when a peak mark held focus (WCAG 2.4.3)", () => {
    // FO-FOCUSRETURN: disarming strips every mark's tabIndex/role, so an Escape
    // exit while a mark holds keyboard focus would drop focus to <body>. The
    // handler re-anchors to the "+ Peak" button -- the keyboard user's stable
    // handle -- before the disarm re-render makes the mark inert.
    const { container } = render(
      <TracePlate
        {...base}
        addPeakArmed
        onToggleAddPeak={() => {}}
        interaction={{ onXDomain: () => {}, onAddPeak: () => {}, onClickPeak: () => {} }}
      />,
    );
    const mark = container.querySelector(
      '[data-role="plot-peaks"] [role="button"]',
    ) as HTMLElement;
    expect(mark).toBeTruthy();
    // Resolve the toolbar button via raw DOM: RTL's accessible-name engine trips
    // (reads `.name` of undefined) when it has to compute a role for the SVG
    // <g role="button"> peak marks, so getByRole/getByText cannot be used here.
    const addPeakBtn = [...container.querySelectorAll("button")].find(
      (b) => b.textContent?.trim() === "+ Peak",
    ) as HTMLElement;
    expect(addPeakBtn).toBeTruthy();
    mark.focus();
    expect(document.activeElement).toBe(mark);
    fireEvent.keyDown(document.body, { key: "Escape" });
    expect(document.activeElement).toBe(addPeakBtn);
  });

  it("Escape-disarm re-anchors focus from the add-at-q field to '+ Peak' (WCAG 2.4.3)", () => {
    // The add-at-q field is armed-only and unmounts on disarm, so an Escape exit
    // while it holds focus would also drop to <body>. (No onClickPeak here, so
    // no SVG role=button marks exist to trip the accessible-name engine.)
    render(
      <TracePlate
        {...base}
        addPeakArmed
        onToggleAddPeak={() => {}}
        interaction={{ onXDomain: () => {}, onAddPeak: () => {} }}
      />,
    );
    const input = screen.getByLabelText("q value for new peak");
    input.focus();
    expect(document.activeElement).toBe(input);
    fireEvent.keyDown(document.body, { key: "Escape" });
    expect(document.activeElement).toBe(
      screen.getByRole("button", { name: "+ Peak" }),
    );
  });

  it("Escape-disarm leaves focus alone when no peak mark held it (no yank)", () => {
    // The re-anchor is scoped to the about-to-go-inert marks: pressing Escape
    // with focus already on the "+ Peak" button (or elsewhere off the plot)
    // must not be hijacked.
    render(
      <TracePlate
        {...base}
        addPeakArmed
        onToggleAddPeak={() => {}}
        interaction={{ onXDomain: () => {}, onAddPeak: () => {} }}
      />,
    );
    const addPeakBtn = screen.getByText("+ Peak").closest("button")!;
    addPeakBtn.focus();
    fireEvent.keyDown(document.body, { key: "Escape" });
    expect(document.activeElement).toBe(addPeakBtn);
  });

  it("Escape while disarmed does nothing (no spurious re-arm toggle)", () => {
    const onToggleAddPeak = vi.fn();
    render(
      <TracePlate
        {...base}
        onToggleAddPeak={onToggleAddPeak}
        interaction={{ onXDomain: () => {}, onAddPeak: () => {} }}
      />,
    );
    fireEvent.keyDown(document.body, { key: "Escape" });
    expect(onToggleAddPeak).not.toHaveBeenCalled();
  });

  it("Escape with an open modal dialog does NOT disarm (the dialog owns Escape)", () => {
    const onToggleAddPeak = vi.fn();
    render(
      <>
        <TracePlate
          {...base}
          addPeakArmed
          onToggleAddPeak={onToggleAddPeak}
          interaction={{ onXDomain: () => {}, onAddPeak: () => {} }}
        />
        {/* What ModalShell emits while open — e.g. the custom-index modal. */}
        <div role="dialog" aria-modal="true" aria-label="stub dialog" />
      </>,
    );
    fireEvent.keyDown(document.body, { key: "Escape" });
    expect(onToggleAddPeak).not.toHaveBeenCalled();
  });

  it("an already-consumed Escape does NOT disarm even with the dialog gone (real-browser closing press)", () => {
    const onToggleAddPeak = vi.fn();
    render(
      <TracePlate
        {...base}
        addPeakArmed
        onToggleAddPeak={onToggleAddPeak}
        interaction={{ onXDomain: () => {}, onAddPeak: () => {} }}
      />,
    );
    // In a real browser a microtask checkpoint runs between ModalShell's
    // document-level close listener and TracePlate's window-level one, so the
    // dialog is already unmounted when the disarm guard fires — the only trace
    // of the modal is the preventDefault stamp ModalShell left on the event.
    const consume = (e: KeyboardEvent): void => e.preventDefault();
    document.addEventListener("keydown", consume);
    try {
      fireEvent.keyDown(document.body, { key: "Escape" });
    } finally {
      document.removeEventListener("keydown", consume);
    }
    expect(onToggleAddPeak).not.toHaveBeenCalled();
  });

  it("Escape inside the add-at-q input disarms (the input has no local Escape behavior)", () => {
    const onToggleAddPeak = vi.fn();
    render(
      <TracePlate
        {...base}
        addPeakArmed
        onToggleAddPeak={onToggleAddPeak}
        interaction={{ onXDomain: () => {}, onAddPeak: () => {} }}
      />,
    );
    const input = screen.getByLabelText("q value for new peak");
    fireEvent.keyDown(input, { key: "Escape" });
    expect(onToggleAddPeak).toHaveBeenCalledTimes(1);
  });

  it("armed hint names the Esc exit when the toggle is wired", () => {
    render(
      <TracePlate
        {...base}
        addPeakArmed
        onToggleAddPeak={() => {}}
        interaction={{ onXDomain: () => {}, onAddPeak: () => {}, onClickPeak: () => {} }}
      />,
    );
    const hint = screen.getByTestId("peak-edit-hint");
    expect(hint.textContent).toContain("Esc exits.");
  });

  it("armed hint omits the Esc sentence when no onToggleAddPeak is wired (controls don't lie)", () => {
    render(
      <TracePlate
        {...base}
        addPeakArmed
        interaction={{ onXDomain: () => {}, onAddPeak: () => {}, onClickPeak: () => {} }}
      />,
    );
    const hint = screen.getByTestId("peak-edit-hint");
    expect(hint.textContent).not.toContain("Esc exits.");
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
