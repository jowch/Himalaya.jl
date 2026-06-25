import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { SpecCell } from "../../src/print/components/SpecCell";
import { KeptCell } from "../../src/print/components/KeptCell";
import { StatusCell } from "../../src/print/components/StatusCell";

// ---------------------------------------------------------------------------
// SpecCell
// ---------------------------------------------------------------------------

describe("<SpecCell> name and id", () => {
  it("renders the sample name", () => {
    render(<SpecCell name="Sample A" sampleId="s-001" />);
    expect(screen.getByText("Sample A")).toBeInTheDocument();
  });

  it("renders the sample id", () => {
    render(<SpecCell name="Sample A" sampleId="s-001" />);
    expect(screen.getByText("s-001")).toBeInTheDocument();
  });

  it("name is in data-role=spec-name", () => {
    const { container } = render(<SpecCell name="Sample A" sampleId="s-001" />);
    expect(container.querySelector("[data-role='spec-name']")).toHaveTextContent("Sample A");
  });

  it("id is in data-role=spec-id", () => {
    const { container } = render(<SpecCell name="Sample A" sampleId="s-001" />);
    expect(container.querySelector("[data-role='spec-id']")).toHaveTextContent("s-001");
  });
});

describe("<SpecCell> slot chip", () => {
  it("renders a slot-chip when slotIndex is provided", () => {
    render(<SpecCell name="Sample A" sampleId="s-001" slotIndex={3} />);
    const chip = screen.getByTestId("slot-chip");
    expect(chip).toBeInTheDocument();
    expect(chip).toHaveTextContent("slot 3");
  });

  it("renders no slot-chip when slotIndex is absent", () => {
    render(<SpecCell name="Sample A" sampleId="s-001" />);
    expect(screen.queryByTestId("slot-chip")).toBeNull();
  });

  it("renders no slot-chip when slotIndex is absent (even when screened)", () => {
    render(<SpecCell name="Sample A" sampleId="s-001" screened />);
    expect(screen.queryByTestId("slot-chip")).toBeNull();
  });

  it("renders the root data-testid=spec-cell", () => {
    render(<SpecCell name="Sample A" sampleId="s-001" />);
    expect(screen.getByTestId("spec-cell")).toBeInTheDocument();
  });
});

// ---------------------------------------------------------------------------
// KeptCell
// ---------------------------------------------------------------------------

describe("<KeptCell> counts", () => {
  it("renders the kept count text", () => {
    render(<KeptCell kept={3} total={5} />);
    expect(screen.getByText("3")).toBeInTheDocument();
  });

  it("renders the total with slash prefix", () => {
    render(<KeptCell kept={3} total={5} />);
    expect(screen.getByText("/ 5")).toBeInTheDocument();
  });

  it("kept count is in data-role=kept-count", () => {
    const { container } = render(<KeptCell kept={3} total={5} />);
    expect(container.querySelector("[data-role='kept-count']")).toHaveTextContent("3");
  });

  it("total is in data-role=kept-total", () => {
    const { container } = render(<KeptCell kept={3} total={5} />);
    expect(container.querySelector("[data-role='kept-total']")).toHaveTextContent("/ 5");
  });

  it("renders the root data-testid=kept-cell", () => {
    render(<KeptCell kept={3} total={5} />);
    expect(screen.getByTestId("kept-cell")).toBeInTheDocument();
  });
});

describe("<KeptCell> dropped span", () => {
  it("renders data-role=kept-dropped with correct text when dropped > 0", () => {
    const { container } = render(<KeptCell kept={4} total={5} dropped={1} />);
    const el = container.querySelector("[data-role='kept-dropped']");
    expect(el).toBeInTheDocument();
    expect(el).toHaveTextContent("1 dropped");
  });

  it("does NOT render data-role=kept-dropped when dropped=0", () => {
    const { container } = render(<KeptCell kept={5} total={5} dropped={0} />);
    expect(container.querySelector("[data-role='kept-dropped']")).not.toBeInTheDocument();
  });

  it("does NOT render data-role=kept-dropped when dropped is omitted", () => {
    const { container } = render(<KeptCell kept={5} total={5} />);
    expect(container.querySelector("[data-role='kept-dropped']")).not.toBeInTheDocument();
  });
});

// ---------------------------------------------------------------------------
// StatusCell
// ---------------------------------------------------------------------------

describe("<StatusCell> with phase", () => {
  it("renders PhaseChip with the phase text when phase is set", () => {
    render(<StatusCell phase="Pn3m" />);
    expect(screen.getByTestId("phase-chip")).toBeInTheDocument();
    expect(screen.getByText("Pn3m")).toBeInTheDocument();
  });

  it("does NOT render data-role=status-unset when phase is provided", () => {
    const { container } = render(<StatusCell phase="Pn3m" />);
    expect(container.querySelector("[data-role='status-unset']")).not.toBeInTheDocument();
  });

  it("renders Im3m phase correctly", () => {
    render(<StatusCell phase="Im3m" />);
    expect(screen.getByText("Im3m")).toBeInTheDocument();
  });

  it("renders the root data-testid=status-cell", () => {
    render(<StatusCell phase="Pn3m" />);
    expect(screen.getByTestId("status-cell")).toBeInTheDocument();
  });
});

describe("<StatusCell> unset states", () => {
  it("renders data-role=status-unset with 'Not indexed' text when phase is null", () => {
    const { container } = render(<StatusCell phase={null} />);
    const el = container.querySelector("[data-role='status-unset']");
    expect(el).toBeInTheDocument();
    expect(el).toHaveTextContent("Not indexed");
  });

  it("renders data-role=status-unset when phase is empty string", () => {
    const { container } = render(<StatusCell phase="" />);
    expect(container.querySelector("[data-role='status-unset']")).toBeInTheDocument();
  });

  it("does NOT render PhaseChip when phase is null", () => {
    render(<StatusCell phase={null} />);
    expect(screen.queryByTestId("phase-chip")).not.toBeInTheDocument();
  });

  it("does NOT render PhaseChip when phase is empty string", () => {
    render(<StatusCell phase="" />);
    expect(screen.queryByTestId("phase-chip")).not.toBeInTheDocument();
  });

  it("renders data-role=status-unset when phase is omitted", () => {
    const { container } = render(<StatusCell />);
    expect(container.querySelector("[data-role='status-unset']")).toBeInTheDocument();
  });

  it("Dot with tone=muted is present in unset state", () => {
    const { container } = render(<StatusCell phase={null} />);
    expect(container.querySelector("[data-tone='muted']")).toBeInTheDocument();
  });
});

describe("<StatusCell> noExposures (SA-ZEROEXP)", () => {
  it("reads 'No exposures' (not 'Not indexed') and carries its own data-role", () => {
    const { container } = render(<StatusCell phase={null} noExposures />);
    const el = container.querySelector("[data-role='status-no-exposures']");
    expect(el).toBeInTheDocument();
    expect(el).toHaveTextContent("No exposures");
    expect(container.querySelector("[data-role='status-unset']")).not.toBeInTheDocument();
  });

  it("takes precedence over a phase (an empty sample can't carry a phase)", () => {
    render(<StatusCell phase="Pn3m" noExposures />);
    expect(screen.getByText("No exposures")).toBeInTheDocument();
    expect(screen.queryByTestId("phase-chip")).not.toBeInTheDocument();
  });
});

describe("<StatusCell> form factor (representative exposure declared form_factor)", () => {
  it("renders a distinct 'Form factor' status, not 'Not indexed'", () => {
    const { container } = render(<StatusCell phase={null} formFactor />);
    const el = container.querySelector("[data-role='status-form-factor']");
    expect(el).toBeInTheDocument();
    expect(el).toHaveTextContent("Form factor");
    expect(container.querySelector("[data-role='status-unset']")).not.toBeInTheDocument();
  });

  it("takes precedence over the door invite (a classified sample is not an 'Index' invitation)", () => {
    const { container } = render(<StatusCell phase={null} formFactor door />);
    expect(container.querySelector("[data-role='status-form-factor']")).toBeInTheDocument();
    expect(container.querySelector("[data-role='status-index-invite']")).not.toBeInTheDocument();
  });

  it("yields to noExposures (an empty sample can't be form factor)", () => {
    const { container } = render(<StatusCell phase={null} formFactor noExposures />);
    expect(container.querySelector("[data-role='status-no-exposures']")).toBeInTheDocument();
    expect(container.querySelector("[data-role='status-form-factor']")).not.toBeInTheDocument();
  });
});

describe("<StatusCell> door invite (SA-RESCORE3 F9)", () => {
  it("an unindexed DOOR reads as an 'Index' invitation, not a dead 'Not indexed' status", () => {
    const { container } = render(<StatusCell phase={null} door />);
    const invite = container.querySelector("[data-role='status-index-invite']");
    expect(invite).toBeInTheDocument();
    expect(invite).toHaveTextContent("Index");
    // The passive status must not also render — the door affords its action.
    expect(container.querySelector("[data-role='status-unset']")).not.toBeInTheDocument();
  });

  it("a non-door unindexed cell keeps the passive 'Not indexed' status (no false invite)", () => {
    const { container } = render(<StatusCell phase={null} />);
    expect(container.querySelector("[data-role='status-unset']")).toBeInTheDocument();
    expect(container.querySelector("[data-role='status-index-invite']")).not.toBeInTheDocument();
  });

  it("an indexed door keeps the PhaseChip (already an obvious open affordance)", () => {
    render(<StatusCell phase="Pn3m" door />);
    expect(screen.getByTestId("phase-chip")).toBeInTheDocument();
    expect(screen.queryByText(/^Index$/)).not.toBeInTheDocument();
  });
});
