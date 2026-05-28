import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { ScopingPhaseStrip } from "../src/components/ScopingPhaseStrip";

describe("ScopingPhaseStrip", () => {
  it("renders one segment per member", () => {
    render(<ScopingPhaseStrip reads={[
      { dominant: "Pn3m", coexist: null },
      { dominant: "Lamellar", coexist: null },
    ]} />);
    expect(screen.getAllByTestId("scoping-ps-seg")).toHaveLength(2);
  });
  it("captions a transition first → last", () => {
    render(<ScopingPhaseStrip reads={[
      { dominant: "Pn3m", coexist: null },
      { dominant: "Lamellar", coexist: null },
    ]} />);
    const cap = screen.getByTestId("scoping-ps-cap");
    expect(cap).toHaveTextContent("Pn3m");
    expect(cap).toHaveTextContent("Lamellar");
  });
  it("captions 'throughout' when first and last agree", () => {
    render(<ScopingPhaseStrip reads={[
      { dominant: "Pn3m", coexist: null },
      { dominant: "Pn3m", coexist: null },
    ]} />);
    expect(screen.getByTestId("scoping-ps-cap")).toHaveTextContent(/throughout/i);
  });
  it("renders a not-yet-indexed caption when no members are indexed", () => {
    render(<ScopingPhaseStrip reads={[{ dominant: null, coexist: null }]} />);
    expect(screen.getByTestId("scoping-ps-cap")).toHaveTextContent(/not yet indexed/i);
  });
});
