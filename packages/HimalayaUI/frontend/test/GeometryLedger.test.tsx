import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { GeometryLedger, type GeometryRow } from "../src/print/components/GeometryLedger";

const ROWS: GeometryRow[] = [
  { key: "beam_energy", label: "Beam energy", value: "9.00 keV", source: "prp" },
  { key: "beam_center_x", label: "Beam center X", value: "421.4 px", source: "setup" },
  { key: "flight_path", label: "Flight path", value: "1.81 m", source: "user" },
];

const cb = { onOverride: () => {}, onRevert: () => {}, onUndo: () => {}, canUndo: false };

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
  it("Override calls onOverride with the key", () => {
    const onOverride = vi.fn();
    render(<GeometryLedger rows={ROWS} {...cb} onOverride={onOverride} />);
    fireEvent.click(screen.getAllByRole("button", { name: /override/i })[0]!);
    expect(onOverride).toHaveBeenCalledWith("beam_energy");
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

  // --- inline edit-in-place tests ---

  it("when editingKey matches a row, that row shows an input seeded with editDraft", () => {
    render(
      <GeometryLedger
        rows={ROWS}
        {...cb}
        editingKey="beam_energy"
        editDraft="9.00 keV"
        onEditDraftChange={() => {}}
        onEditCommit={() => {}}
        onEditCancel={() => {}}
      />
    );
    const inp = screen.getByRole("textbox", { name: /override beam energy/i });
    expect(inp).toBeInTheDocument();
    expect((inp as HTMLInputElement).value).toBe("9.00 keV");
  });

  it("onEditDraftChange fires when the inline input changes", () => {
    const onChange = vi.fn();
    render(
      <GeometryLedger
        rows={ROWS}
        {...cb}
        editingKey="beam_energy"
        editDraft="9.00 keV"
        onEditDraftChange={onChange}
        onEditCommit={() => {}}
        onEditCancel={() => {}}
      />
    );
    const inp = screen.getByRole("textbox", { name: /override beam energy/i });
    fireEvent.change(inp, { target: { value: "10.00 keV" } });
    expect(onChange).toHaveBeenCalledWith("10.00 keV");
  });

  it("onEditCommit fires on Enter in the inline input", () => {
    const onCommit = vi.fn();
    render(
      <GeometryLedger
        rows={ROWS}
        {...cb}
        editingKey="beam_energy"
        editDraft="9.00 keV"
        onEditDraftChange={() => {}}
        onEditCommit={onCommit}
        onEditCancel={() => {}}
      />
    );
    const inp = screen.getByRole("textbox", { name: /override beam energy/i });
    fireEvent.keyDown(inp, { key: "Enter" });
    expect(onCommit).toHaveBeenCalled();
  });

  it("onEditCommit fires on blur from the inline input", () => {
    const onCommit = vi.fn();
    render(
      <GeometryLedger
        rows={ROWS}
        {...cb}
        editingKey="beam_energy"
        editDraft="9.00 keV"
        onEditDraftChange={() => {}}
        onEditCommit={onCommit}
        onEditCancel={() => {}}
      />
    );
    const inp = screen.getByRole("textbox", { name: /override beam energy/i });
    fireEvent.blur(inp);
    expect(onCommit).toHaveBeenCalled();
  });

  it("onEditCancel fires on Escape in the inline input", () => {
    const onCancel = vi.fn();
    render(
      <GeometryLedger
        rows={ROWS}
        {...cb}
        editingKey="beam_energy"
        editDraft="9.00 keV"
        onEditDraftChange={() => {}}
        onEditCommit={() => {}}
        onEditCancel={onCancel}
      />
    );
    const inp = screen.getByRole("textbox", { name: /override beam energy/i });
    fireEvent.keyDown(inp, { key: "Escape" });
    expect(onCancel).toHaveBeenCalled();
  });

  it("non-editing rows still show their text value when a different row is being edited", () => {
    render(
      <GeometryLedger
        rows={ROWS}
        {...cb}
        editingKey="beam_energy"
        editDraft="9.00 keV"
        onEditDraftChange={() => {}}
        onEditCommit={() => {}}
        onEditCancel={() => {}}
      />
    );
    expect(screen.getByText("421.4 px")).toBeInTheDocument();
    expect(screen.getByText("1.81 m")).toBeInTheDocument();
  });
});
