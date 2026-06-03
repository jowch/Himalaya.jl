import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { ScopePlate } from "../../src/print/components/ScopePlate";
import type { PhaseSegment } from "../../src/print/ui";

const PREVIEW: PhaseSegment[] = [{ phase: "Pn3m" }, { phase: "Lamellar" }];
const base = {
  seriesName: "LL37 titration of lipid 1-2",
  grouping: <>You selected <strong>6 samples</strong>.</>,
  orderedBy: "LL37 : lipid ratio",
  orderNote: "Read from the sample names.",
  count: "6 samples · low to high",
  rows: <div data-testid="rows-slot" />,
  candidates: <div data-testid="cands-slot" />,
  preview: PREVIEW,
  footNote: "Confirming records the ratio on every sample.",
} as const;

describe("<ScopePlate>", () => {
  it("renders the title, kicker, ordered-by value + note, count, and both slots", () => {
    render(<ScopePlate {...base} footState={{ kind: "ready", text: "All 6 values confirmed — ready to build" }} />);
    expect(screen.getByText("LL37 titration of lipid 1-2")).toBeInTheDocument();
    expect(screen.getByText(/New series/i)).toBeInTheDocument();
    expect(screen.getByTestId("order-field")).toHaveTextContent("LL37 : lipid ratio");
    expect(screen.getByText("Read from the sample names.")).toBeInTheDocument();
    expect(screen.getByText("6 samples · low to high")).toBeInTheDocument();
    expect(screen.getByTestId("rows-slot")).toBeInTheDocument();
    expect(screen.getByTestId("cands-slot")).toBeInTheDocument();
    expect(screen.getAllByTestId("ps-seg").length).toBe(2);
  });
  it("gates the build button when buildDisabled", () => {
    render(<ScopePlate {...base} buildDisabled footState={{ kind: "warn", text: "1 value to check before you can build" }} />);
    expect(screen.getByText(/1 value to check/)).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /confirm & build/i })).toBeDisabled();
  });
  it("fires onBuild when enabled, and shows Undo only when onUndo is given", () => {
    const onBuild = vi.fn();
    const { rerender } = render(<ScopePlate {...base} footState={{ kind: "ready", text: "ready" }} onBuild={onBuild} />);
    expect(screen.queryByRole("button", { name: /undo/i })).toBeNull();
    fireEvent.click(screen.getByRole("button", { name: /confirm & build/i }));
    expect(onBuild).toHaveBeenCalledOnce();
    rerender(<ScopePlate {...base} footState={{ kind: "ready", text: "ready" }} onUndo={() => {}} undoLabel="resolved smp_04" />);
    expect(screen.getByRole("button", { name: /undo/i })).toBeInTheDocument();
  });
});
