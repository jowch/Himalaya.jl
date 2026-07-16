import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { CandidateRow, CandidateList } from "../../src/print/components/CandidateRow";

describe("<CandidateRow> rendering", () => {
  it("renders the phase name via PhaseChip text", () => {
    render(<CandidateRow phase="Pn3m" score={0.91} why="dominant — three sharp reflections" />);
    expect(screen.getByText("Pn3m")).toBeInTheDocument();
  });

  it("renders the why text", () => {
    render(<CandidateRow phase="Pn3m" score={0.91} why="dominant — three sharp reflections" />);
    expect(screen.getByText("dominant — three sharp reflections")).toBeInTheDocument();
  });

  it("renders the numeric score formatted to 2 decimal places", () => {
    render(<CandidateRow phase="Pn3m" score={0.91} why="three peaks" />);
    expect(screen.getByText("0.91")).toBeInTheDocument();
  });

  it('renders "custom" when score is null', () => {
    render(<CandidateRow phase="Pn3m" score={null} why="hand-drawn index" />);
    expect(screen.getByText("custom")).toBeInTheDocument();
  });

  it("renders a ScoreBar ([data-score-bar]) when score is not null", () => {
    const { container } = render(<CandidateRow phase="Pn3m" score={0.91} why="three peaks" />);
    expect(container.querySelector("[data-score-bar]")).toBeInTheDocument();
  });

  it("does NOT render a ScoreBar when score is null", () => {
    const { container } = render(<CandidateRow phase="Pn3m" score={null} why="hand-drawn" />);
    expect(container.querySelector("[data-score-bar]")).not.toBeInTheDocument();
  });
});

describe("<CandidateRow> BonnetBadge", () => {
  it("renders BonnetBadge when bonnet=true", () => {
    render(<CandidateRow phase="Im3m" score={0.84} why="explains two subtle peaks" bonnet />);
    expect(screen.getByTestId("bonnet-badge")).toBeInTheDocument();
  });

  it("does NOT render BonnetBadge when bonnet is absent", () => {
    render(<CandidateRow phase="Im3m" score={0.84} why="three peaks" />);
    expect(screen.queryByTestId("bonnet-badge")).not.toBeInTheDocument();
  });
});

describe("<CandidateRow> interactivity", () => {
  it("is rendered as a button", () => {
    render(<CandidateRow phase="Pn3m" score={0.91} why="three peaks" />);
    expect(screen.getByRole("button")).toBeInTheDocument();
  });

  it("calls onToggle when clicked", () => {
    const onToggle = vi.fn();
    render(<CandidateRow phase="Pn3m" score={0.91} why="three peaks" onToggle={onToggle} />);
    fireEvent.click(screen.getByRole("button"));
    expect(onToggle).toHaveBeenCalledTimes(1);
  });

  it("does not throw when clicked without onToggle", () => {
    render(<CandidateRow phase="Pn3m" score={0.91} why="three peaks" />);
    expect(() => fireEvent.click(screen.getByRole("button"))).not.toThrow();
  });
});

describe("<CandidateRow> selected state", () => {
  it("sets aria-pressed=true when selected=true", () => {
    render(<CandidateRow phase="Pn3m" score={0.91} why="three peaks" selected />);
    expect(screen.getByRole("button")).toHaveAttribute("aria-pressed", "true");
  });

  it("sets aria-pressed=false when selected=false", () => {
    render(<CandidateRow phase="Pn3m" score={0.91} why="three peaks" selected={false} />);
    expect(screen.getByRole("button")).toHaveAttribute("aria-pressed", "false");
  });

  it("sets aria-pressed=false when selected is absent", () => {
    render(<CandidateRow phase="Pn3m" score={0.91} why="three peaks" />);
    expect(screen.getByRole("button")).toHaveAttribute("aria-pressed", "false");
  });

  it("passes data-selected to the Card when selected=true", () => {
    render(<CandidateRow phase="Pn3m" score={0.91} why="three peaks" selected />);
    expect(screen.getByRole("button")).toHaveAttribute("data-selected", "true");
  });

  it("CheckCircle reflects selected=true via data-checked", () => {
    render(<CandidateRow phase="Pn3m" score={0.91} why="three peaks" selected />);
    expect(screen.getByTestId("check-circle")).toHaveAttribute("data-checked", "true");
  });

  it("CheckCircle reflects selected=false (no data-checked)", () => {
    render(<CandidateRow phase="Pn3m" score={0.91} why="three peaks" />);
    // data-checked is set to undefined (not present) when unchecked
    expect(screen.getByTestId("check-circle")).not.toHaveAttribute("data-checked");
  });
});

describe("<CandidateList>", () => {
  it("renders its CandidateRow children", () => {
    render(
      <CandidateList>
        <CandidateRow phase="Pn3m" score={0.91} why="three peaks" />
        <CandidateRow phase="Im3m" score={0.84} why="explains two" />
      </CandidateList>
    );
    expect(screen.getByText("Pn3m")).toBeInTheDocument();
    expect(screen.getByText("Im3m")).toBeInTheDocument();
  });
});
