import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, fireEvent, render } from "@testing-library/react";
import { AssignmentCart } from "../src/components/AssignmentCart";
import type { Assignment, IndexEntry } from "../src/api";

function ix(id: number, phase: string, opts: Partial<IndexEntry> = {}): IndexEntry {
  return {
    id, exposure_id: 1, phase, basis: 0.1, score: 0.9, r_squared: 0.99,
    lattice_d: 197, ngc: null, status: "candidate", kind: "auto",
    inputs_hash: null,
    peaks: [{ peak_id: id, ratio_position: 1, residual: 0, q_observed: 0.045 }],
    predicted_q: [0.045], ...opts,
  };
}

const INDICES = [ix(10, "Pn3m"), ix(11, "Im3m", { lattice_d: 252 })];

interface Harness {
  onAdd: ReturnType<typeof vi.fn>;
  onRemove: ReturnType<typeof vi.fn>;
  onSetState: ReturnType<typeof vi.fn>;
  onOpenCustom: ReturnType<typeof vi.fn>;
}

function renderCart(assignment: Assignment, peakCount = 1): Harness {
  const onAdd = vi.fn(), onRemove = vi.fn(), onSetState = vi.fn(), onOpenCustom = vi.fn();
  render(
    <AssignmentCart
      assignment={assignment}
      indices={INDICES}
      peakCount={peakCount}
      latticeUnit="Å"
      onRemovePhase={onRemove}
      onSetState={onSetState}
      onOpenCustom={onOpenCustom}
    />,
  );
  return { onAdd, onRemove, onSetState, onOpenCustom };
}

describe("<AssignmentCart> (Plan D-4)", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("renders a phase block per member for the indexed state", () => {
    renderCart({ exposure_id: 1, state: "indexed", members: [10, 11] });
    expect(screen.getByTestId("assignment-block-10")).toBeInTheDocument();
    expect(screen.getByTestId("assignment-block-11")).toBeInTheDocument();
  });

  it("× removes a phase via onRemovePhase", () => {
    const h = renderCart({ exposure_id: 1, state: "indexed", members: [10] });
    fireEvent.click(screen.getByLabelText(/remove Pn3m/i));
    expect(h.onRemove).toHaveBeenCalledWith(10);
  });

  it("renders the empty/form-factor declaration when indexed with 0 members", () => {
    renderCart({ exposure_id: 1, state: "indexed", members: [] });
    expect(screen.getByTestId("assignment-empty")).toBeInTheDocument();
  });

  it("renders the form-factor declaration chip for form_factor state", () => {
    renderCart({ exposure_id: 1, state: "form_factor", members: [] });
    expect(screen.getByTestId("assignment-declaration")).toHaveTextContent(/form factor/i);
  });

  it("renders the no-scattering declaration for null state", () => {
    renderCart({ exposure_id: 1, state: "null", members: [] });
    expect(screen.getByTestId("assignment-declaration")).toHaveTextContent(/no scattering/i);
  });

  it("SegmentedControl switches state via onSetState", () => {
    const h = renderCart({ exposure_id: 1, state: "indexed", members: [10] });
    fireEvent.click(screen.getByRole("radio", { name: /form factor/i }));
    expect(h.onSetState).toHaveBeenCalledWith("form_factor");
  });

  it("footer custom-index button opens the modal and reads as a button", () => {
    const h = renderCart({ exposure_id: 1, state: "indexed", members: [10] });
    const foot = screen.getByTestId("custom-index-footer");
    expect(foot.tagName).toBe("BUTTON");
    fireEvent.click(foot);
    expect(h.onOpenCustom).toHaveBeenCalled();
  });

  it("Bonnet note: cubic in cart + no second cubic + unexplained peaks → suggestion", () => {
    // Pn3m claims 1 peak; peakCount=3 → 2 unexplained → suggestion note shows.
    renderCart({ exposure_id: 1, state: "indexed", members: [10] }, 3);
    expect(screen.getByTestId("bonnet-note")).toHaveTextContent(/unexplained|account/i);
  });

  it("Bonnet note: two coexisting cubics → consistency note shows", () => {
    renderCart({ exposure_id: 1, state: "indexed", members: [10, 11] }, 2);
    expect(screen.getByTestId("bonnet-note")).toHaveTextContent(/coexist|ratio|consistent/i);
  });

  it("Bonnet note absent when not substantive (single non-cubic, all peaks explained)", () => {
    renderCart({ exposure_id: 1, state: "indexed", members: [10] }, 1);
    expect(screen.queryByTestId("bonnet-note")).not.toBeInTheDocument();
  });
});
