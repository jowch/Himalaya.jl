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
});
