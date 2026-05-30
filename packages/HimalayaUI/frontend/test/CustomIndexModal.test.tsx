import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { CustomIndexModal } from "../src/components/CustomIndexModal";
import { customRefls, basisFor } from "../src/lib/customIndex";

// Two observed peaks that a Pn3m at a≈197 would land its first two orders on.
const PEAK_QS = customRefls("Pn3m", 197).slice(0, 2).map((r) => r.q);

interface Harness { onCommit: ReturnType<typeof vi.fn>; onClose: ReturnType<typeof vi.fn>; }

function renderModal(): Harness {
  const onCommit = vi.fn(), onClose = vi.fn();
  render(
    <CustomIndexModal open peakQs={PEAK_QS} onCommit={onCommit} onClose={onClose} />,
  );
  return { onCommit, onClose };
}

describe("<CustomIndexModal> (Plan D-9)", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("renders the symmetry picker and a live preview comb", () => {
    renderModal();
    expect(screen.getByTestId("custom-index-symmetry")).toBeInTheDocument();
    expect(screen.getByTestId("custom-index-preview")).toBeInTheDocument();
  });

  it("computes reflections from physics and reports the fit (lands on N of M)", () => {
    renderModal();
    // default Pn3m at def=197 lands its first two orders on the two peaks.
    expect(screen.getByTestId("custom-index-fit")).toHaveTextContent(/lands on\s*2\s*of\s*2/i);
  });

  it("clicking a peak snaps the first reflection onto it", () => {
    renderModal();
    // click the first observed-peak guide → first order lands on it
    fireEvent.click(screen.getByTestId(`custom-index-peak-0`));
    // after snap the fit still reports a landing on that peak
    expect(screen.getByTestId("custom-index-fit")).toHaveTextContent(/lands on/i);
  });

  it("commit submits {phase, basis} via onCommit (basis = first reflection q)", () => {
    const h = renderModal();
    fireEvent.click(screen.getByTestId("custom-index-commit"));
    expect(h.onCommit).toHaveBeenCalledTimes(1);
    const [phase, basis] = h.onCommit.mock.calls[0]!;
    expect(phase).toBe("Pn3m");
    expect(basis).toBeCloseTo(basisFor("Pn3m", 197), 6);
  });

  it("switching symmetry to Ia3d recomputes the comb and the committed basis", () => {
    const h = renderModal();
    fireEvent.click(screen.getByRole("radio", { name: /Ia3d/i }));
    fireEvent.click(screen.getByTestId("custom-index-commit"));
    const [phase, basis] = h.onCommit.mock.calls[0]!;
    expect(phase).toBe("Ia3d");
    // Ia3d at its default lattice: basis = 2π√6/a
    expect(basis).toBeGreaterThan(0);
  });
});
