import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { GeometryLedger, type GeometryRow } from "../src/print/components/GeometryLedger";

const ROWS: GeometryRow[] = [
  { key: "beam_energy", label: "Beam energy", value: "9.00 keV", source: "prp", editValue: "9.00" },
  { key: "beam_center_x", label: "Beam center X", value: "421.4 px", source: "setup", editValue: "421.4" },
  { key: "flight_path", label: "Flight path", value: "1.81 m", source: "user", editValue: "1.81" },
];

const cb = {
  onCommit: () => {},
  onRevert: () => {},
  onUndo: () => {},
  canUndo: false,
  onRedo: () => {},
  canRedo: false,
};

describe("GeometryLedger", () => {
  it("renders each row with a provenance chip", () => {
    render(<GeometryLedger rows={ROWS} {...cb} />);
    expect(screen.getByText("Beam energy")).toBeInTheDocument();
    expect(screen.getByText("PRP")).toBeInTheDocument();
    expect(screen.getByText("setup files")).toBeInTheDocument();
    expect(screen.getByText("edited")).toBeInTheDocument();
  });
  it("user rows expose Revert; calls onRevert with the key", () => {
    const onRevert = vi.fn();
    render(<GeometryLedger rows={ROWS} {...cb} onRevert={onRevert} />);
    fireEvent.click(screen.getByRole("button", { name: /revert/i }));
    expect(onRevert).toHaveBeenCalledWith("flight_path");
  });
  it("shows the discrepancy banner when discrepancies > 0", () => {
    render(<GeometryLedger rows={ROWS} {...cb} discrepancyCount={2} />);
    expect(screen.getByText(/geometry check found 2 issues/i)).toBeInTheDocument();
  });
  it("Undo last change is gated on canUndo", () => {
    const onUndo = vi.fn();
    const { rerender } = render(<GeometryLedger rows={ROWS} {...cb} canUndo={false} onUndo={onUndo} />);
    expect(screen.queryByRole("button", { name: /undo last change/i })).toBeNull();
    rerender(<GeometryLedger rows={ROWS} {...cb} canUndo onUndo={onUndo} />);
    fireEvent.click(screen.getByRole("button", { name: /undo last change/i }));
    expect(onUndo).toHaveBeenCalled();
  });
  it("Redo is gated on canRedo", () => {
    const onRedo = vi.fn();
    const { rerender } = render(<GeometryLedger rows={ROWS} {...cb} canRedo={false} onRedo={onRedo} />);
    expect(screen.queryByRole("button", { name: /redo last change/i })).toBeNull();
    rerender(<GeometryLedger rows={ROWS} {...cb} canRedo onRedo={onRedo} />);
    fireEvent.click(screen.getByRole("button", { name: /redo last change/i }));
    expect(onRedo).toHaveBeenCalled();
  });

  // --- inline edit-in-place tests (GeometryLedger owns editing state) ---

  it("Override opens an input on that row, seeded with the row's editValue", () => {
    render(<GeometryLedger rows={ROWS} {...cb} />);
    fireEvent.click(screen.getByRole("button", { name: /override beam energy/i }));
    const inp = screen.getByRole("textbox", { name: /override beam energy/i });
    expect((inp as HTMLInputElement).value).toBe("9.00");
  });

  it("a row with no editValue is not editable -- Override is a no-op", () => {
    const rows: GeometryRow[] = [{ key: "q_units", label: "q units", value: "—", source: "default" }];
    render(<GeometryLedger rows={rows} {...cb} />);
    fireEvent.click(screen.getByRole("button", { name: /override q units/i }));
    expect(screen.queryByRole("textbox", { name: /override q units/i })).toBeNull();
  });

  it("Enter commits the current draft via onCommit(key, value)", () => {
    const onCommit = vi.fn();
    render(<GeometryLedger rows={ROWS} {...cb} onCommit={onCommit} />);
    fireEvent.click(screen.getByRole("button", { name: /override beam energy/i }));
    const inp = screen.getByRole("textbox", { name: /override beam energy/i });
    fireEvent.change(inp, { target: { value: "10.00" } });
    fireEvent.keyDown(inp, { key: "Enter" });
    expect(onCommit).toHaveBeenCalledWith("beam_energy", "10.00");
  });

  it("blur commits via onCommit", () => {
    const onCommit = vi.fn();
    render(<GeometryLedger rows={ROWS} {...cb} onCommit={onCommit} />);
    fireEvent.click(screen.getByRole("button", { name: /override beam energy/i }));
    fireEvent.blur(screen.getByRole("textbox", { name: /override beam energy/i }));
    expect(onCommit).toHaveBeenCalledWith("beam_energy", "9.00");
  });

  it("Enter then a trailing blur commits exactly once (no double undo push)", () => {
    const onCommit = vi.fn();
    render(<GeometryLedger rows={ROWS} {...cb} onCommit={onCommit} />);
    fireEvent.click(screen.getByRole("button", { name: /override beam energy/i }));
    const inp = screen.getByRole("textbox", { name: /override beam energy/i });
    fireEvent.keyDown(inp, { key: "Enter" });
    fireEvent.blur(inp);
    expect(onCommit).toHaveBeenCalledTimes(1);
  });

  it("switching rows mid-edit commits the open row exactly once, then opens the next", () => {
    // Real browsers blur the open input when focus leaves to click another row's
    // Override. The ref guard must serialize: blur A commits once, opening B
    // must not re-fire A's commit (which would risk a stray second undo entry).
    const onCommit = vi.fn();
    render(<GeometryLedger rows={ROWS} {...cb} onCommit={onCommit} />);
    fireEvent.click(screen.getByRole("button", { name: /override beam energy/i }));
    fireEvent.blur(screen.getByRole("textbox", { name: /override beam energy/i }));
    fireEvent.click(screen.getByRole("button", { name: /override beam center x/i }));
    expect(screen.getByRole("textbox", { name: /override beam center x/i })).toBeInTheDocument();
    expect(onCommit).toHaveBeenCalledTimes(1);
    expect(onCommit).toHaveBeenCalledWith("beam_energy", "9.00");
  });

  it("Escape exits edit mode without committing", () => {
    const onCommit = vi.fn();
    render(<GeometryLedger rows={ROWS} {...cb} onCommit={onCommit} />);
    fireEvent.click(screen.getByRole("button", { name: /override beam energy/i }));
    const inp = screen.getByRole("textbox", { name: /override beam energy/i });
    fireEvent.keyDown(inp, { key: "Escape" });
    expect(onCommit).not.toHaveBeenCalled();
    expect(screen.queryByRole("textbox", { name: /override beam energy/i })).toBeNull();
    expect(screen.getByText("9.00 keV")).toBeInTheDocument();
  });

  it("source='computed' renders a non-empty provenance chip label", () => {
    const computedRows: GeometryRow[] = [
      { key: "beam_center_x", label: "Beam center X", value: "421.4 px", source: "computed" },
    ];
    render(<GeometryLedger rows={computedRows} {...cb} />);
    expect(screen.getByText("computed")).toBeInTheDocument();
  });

  it("non-editing rows still show their text value when a different row is being edited", () => {
    render(<GeometryLedger rows={ROWS} {...cb} />);
    fireEvent.click(screen.getByRole("button", { name: /override beam energy/i }));
    expect(screen.getByText("421.4 px")).toBeInTheDocument();
    expect(screen.getByText("1.81 m")).toBeInTheDocument();
  });
});
