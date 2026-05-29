import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { PhaseStrip } from "../src/components/ui/PhaseStrip";
import type { PhaseRead } from "../src/lib/scoping/dominantPhase";

// Mirrors the SeriesScopingPage call site: PhaseRead[] → PhaseSegment[].
function strip(reads: PhaseRead[]) {
  return render(
    <PhaseStrip
      size="sm"
      emptyLabel="Members not yet indexed; phase preview unavailable."
      segments={reads.map((r) => ({ phase: r.dominant, coexistWith: r.coexist }))}
    />,
  );
}

describe("Scoping preview strip (PhaseStrip, size=sm)", () => {
  it("renders one segment per member", () => {
    strip([
      { dominant: "Pn3m", coexist: null },
      { dominant: "Lamellar", coexist: null },
    ]);
    expect(screen.getAllByTestId("ps-seg")).toHaveLength(2);
  });

  it("captions a transition first → last", () => {
    strip([
      { dominant: "Pn3m", coexist: null },
      { dominant: "Lamellar", coexist: null },
    ]);
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent("Pn3m");
    expect(cap).toHaveTextContent("Lamellar");
  });

  it("captions 'throughout' when every indexed segment is one phase", () => {
    strip([
      { dominant: "Pn3m", coexist: null },
      { dominant: "Pn3m", coexist: null },
    ]);
    expect(screen.getByTestId("ps-cap")).toHaveTextContent(/throughout/i);
  });

  it("renders the not-yet-indexed empty label when no members are indexed", () => {
    strip([{ dominant: null, coexist: null }]);
    expect(screen.getByTestId("ps-cap")).toHaveTextContent(/not yet indexed/i);
  });
});
