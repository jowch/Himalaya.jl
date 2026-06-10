import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { ScopeCandidateRow } from "../../src/print/components/ScopeCandidateRow";

const TRACE = { q: [0.02, 0.04, 0.08, 0.16], I: [10, 80, 30, 5] };

describe("<ScopeCandidateRow>", () => {
  it("renders the name, the why, and a dimmed sparkline", () => {
    render(<ScopeCandidateRow name="Lipid 1-1 + LL37 1:1"
      why={<>has LL37, but the <strong>1-1 lipid line</strong> — its own series?</>}
      trace={TRACE} phase="Pn3m" />);
    const row = screen.getByTestId("scope-candidate-row");
    expect(row).toHaveTextContent("Lipid 1-1 + LL37 1:1");
    expect(row).toHaveTextContent("1-1 lipid line");
    expect(screen.getByTestId("sparkline")).toBeInTheDocument();
  });
  it("renders the optional smp_ sample identity next to the name", () => {
    render(<ScopeCandidateRow name="Lipid 1-1 + LL37 1:1" sampleId="smp_42"
      why="has LL37" trace={TRACE} />);
    expect(screen.getByText("smp_42")).toBeInTheDocument();
  });
  it("omits the identity line when sampleId is not given", () => {
    render(<ScopeCandidateRow name="x" why="y" trace={TRACE} />);
    expect(screen.queryByText(/^smp_/)).toBeNull();
  });
  it("fires onAdd from the Add-to-series button", () => {
    const onAdd = vi.fn();
    render(<ScopeCandidateRow name="x" why="y" trace={TRACE} onAdd={onAdd} />);
    fireEvent.click(screen.getByRole("button", { name: /add to series/i }));
    expect(onAdd).toHaveBeenCalledOnce();
  });
});
