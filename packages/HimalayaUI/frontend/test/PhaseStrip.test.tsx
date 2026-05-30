import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { PhaseStrip } from "../src/components/ui/PhaseStrip";
import type { PhaseSegment } from "../src/components/ui/PhaseStrip";

const seg = (phase: string | null, coexistWith: string | null = null): PhaseSegment => ({
  phase,
  coexistWith,
});

describe("PhaseStrip", () => {
  it("renders one ps-seg per segment", () => {
    render(<PhaseStrip segments={[seg("Pn3m"), seg("Lamellar"), seg("Im3m")]} />);
    expect(screen.getAllByTestId("ps-seg")).toHaveLength(3);
  });

  it("captions a transition with both phase names and a decorative arrow", () => {
    render(<PhaseStrip segments={[seg("Pn3m"), seg("Lamellar")]} />);
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent("Pn3m");
    expect(cap).toHaveTextContent("Lamellar");
    // arrow glyph present...
    expect(cap).toHaveTextContent("→");
    // ...but decorative (aria-hidden), so not part of the cap's accessible text channel.
    const arrow = cap.querySelector('[aria-hidden="true"]');
    expect(arrow).not.toBeNull();
    expect(arrow).toHaveTextContent("→");
  });

  it("captions a single distinct phase as '<phase> throughout'", () => {
    render(<PhaseStrip segments={[seg("Pn3m"), seg("Pn3m")]} />);
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent("Pn3m");
    expect(cap).toHaveTextContent(/throughout/i);
  });

  it("treats a non-monotone strip as a transition (distinct-count, not first===last)", () => {
    // [Pn3m, Lamellar, Pn3m] has two distinct indexed phases → transition, NOT "throughout".
    render(<PhaseStrip segments={[seg("Pn3m"), seg("Lamellar"), seg("Pn3m")]} />);
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent(/→/);
    expect(cap).not.toHaveTextContent(/throughout/i);
  });

  it("renders the default empty label when no segment is indexed", () => {
    render(<PhaseStrip segments={[seg(null), seg(null)]} />);
    expect(screen.getByTestId("ps-cap")).toHaveTextContent(/no clear phase/i);
  });

  it("honors an emptyLabel override", () => {
    render(
      <PhaseStrip
        segments={[seg(null)]}
        emptyLabel="Members not yet indexed; phase preview unavailable."
      />,
    );
    expect(screen.getByTestId("ps-cap")).toHaveTextContent(/not yet indexed/i);
  });

  it("labels each segment with its phase via aria-label and title", () => {
    render(<PhaseStrip segments={[seg("Pn3m"), seg(null)]} />);
    const segs = screen.getAllByTestId("ps-seg");
    expect(segs[0]).toHaveAttribute("aria-label", "Pn3m");
    expect(segs[0]).toHaveAttribute("title", "Pn3m");
    expect(segs[1]).toHaveAttribute("aria-label", "Unindexed");
    expect(segs[1]).toHaveAttribute("title", "Unindexed");
  });

  it("labels a coexistence segment as 'A + B (coexistence)'", () => {
    render(<PhaseStrip segments={[seg("Pn3m", "Lamellar")]} />);
    const segment = screen.getAllByTestId("ps-seg")[0];
    expect(segment).toHaveAttribute("aria-label", "Pn3m + Lamellar (coexistence)");
    expect(segment).toHaveAttribute("title", "Pn3m + Lamellar (coexistence)");
  });

  it("exposes a data-size attribute reflecting the size prop (default md)", () => {
    const { rerender, container } = render(<PhaseStrip segments={[seg("Pn3m")]} />);
    expect(container.firstChild).toHaveAttribute("data-size", "md");
    rerender(<PhaseStrip segments={[seg("Pn3m")]} size="sm" />);
    expect(container.firstChild).toHaveAttribute("data-size", "sm");
  });

  it("renders a hollow dashed cell for a form_factor segment", () => {
    render(<PhaseStrip segments={[{ phase: null, state: "form_factor" }]} />);
    const segment = screen.getAllByTestId("ps-seg")[0]!;
    expect(segment).toHaveAttribute("data-state", "form_factor");
    expect(segment).toHaveAttribute("aria-label", expect.stringMatching(/form factor/i));
  });

  it("renders a null-state cell distinct from a form-factor cell", () => {
    render(<PhaseStrip segments={[{ phase: null, state: "null" }]} />);
    const segment = screen.getAllByTestId("ps-seg")[0]!;
    expect(segment).toHaveAttribute("data-state", "null");
    // distinct data-state from form_factor
    expect(segment.getAttribute("data-state")).not.toBe("form_factor");
  });

  it("applies the placement className to the root", () => {
    const { container } = render(<PhaseStrip segments={[seg("Pn3m")]} className="mt-5" />);
    expect(container.firstChild).toHaveClass("mt-5");
  });
});
