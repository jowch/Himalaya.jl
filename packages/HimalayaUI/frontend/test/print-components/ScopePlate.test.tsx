import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { ScopePlate } from "../../src/print/components/ScopePlate";
import type { PhaseSegment } from "../../src/print/ui";

const PREVIEW: PhaseSegment[] = [{ phase: "Pn3m" }, { phase: "Lamellar" }];
const base = {
  seriesName: "LL37 titration of lipid 1-2",
  grouping: <>You selected <strong>6 samples</strong>.</>,
  orderedBy: "LL37 : lipid ratio",
  orderNote: "Each value is the sample's stored ratio tag.",
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
    expect(screen.getByText("Each value is the sample's stored ratio tag.")).toBeInTheDocument();
    expect(screen.getByText("6 samples · low to high")).toBeInTheDocument();
    expect(screen.getByTestId("rows-slot")).toBeInTheDocument();
    expect(screen.getByTestId("cands-slot")).toBeInTheDocument();
    expect(screen.getAllByTestId("ps-seg").length).toBe(2);
  });
  it("exposes a sound heading tree with a populated preview: one h1 (series name) + four h2 section labels (the Preview h2 is conditional on segments)", () => {
    render(<ScopePlate {...base} footState={{ kind: "ready", text: "ready" }} />);
    // The series name is the single level-1 heading.
    const h1s = screen.getAllByRole("heading", { level: 1 });
    expect(h1s).toHaveLength(1);
    expect(h1s[0]).toHaveAccessibleName("LL37 titration of lipid 1-2");
    // The four section labels are level-2 headings nested under it.
    expect(screen.getByRole("heading", { level: 2, name: "Ordered by" })).toBeInTheDocument();
    expect(screen.getByRole("heading", { level: 2, name: "The series" })).toBeInTheDocument();
    expect(
      screen.getByRole("heading", { level: 2, name: "Himalaya also found" }),
    ).toBeInTheDocument();
    expect(
      screen.getByRole("heading", { level: 2, name: "Preview · phase across the series" }),
    ).toBeInTheDocument();
  });
  it("names the order-field control with its label so it is not a bare value (WCAG 4.1.2)", () => {
    render(
      <ScopePlate
        {...base}
        orderOptions={[
          { value: "LL37 : lipid ratio", label: "LL37 : lipid ratio" },
          { value: "Time", label: "Time" },
        ]}
        footState={{ kind: "ready", text: "ready" }}
      />,
    );
    // With options the field is a real combobox-style trigger button; its
    // accessible name must read label + value, not just the value.
    expect(
      screen.getByRole("button", { name: /ordered by\s+LL37 : lipid ratio/i }),
    ).toBeInTheDocument();
  });
  it("names the ordering dropdown 'Ordered by' so it is not the generic menu default, with checked option", () => {
    render(
      <ScopePlate
        {...base}
        orderOptions={[
          { value: "LL37 : lipid ratio", label: "LL37 : lipid ratio" },
          { value: "Time", label: "Time" },
        ]}
        onOrderSelect={() => {}}
        footState={{ kind: "ready", text: "ready" }}
      />,
    );
    fireEvent.click(screen.getByTestId("order-field"));
    // The menu's accessible name is the ordering label, not "Choose an option".
    expect(screen.getByRole("menu", { name: "Ordered by" })).toBeInTheDocument();
    // Value-selector → radios; the current value is aria-checked.
    expect(
      screen.getByRole("menuitemradio", { name: "LL37 : lipid ratio" }).getAttribute("aria-checked"),
    ).toBe("true");
    expect(screen.getByRole("menuitemradio", { name: "Time" }).getAttribute("aria-checked")).toBe(
      "false",
    );
  });
  it("omits the preview section entirely when there are no segments (every member skipped)", () => {
    // A zero-segment PhaseStrip would still paint an empty bar + the
    // "No clear phase" caption — a visible artifact previewing nothing. The
    // plate drops the whole section (heading included); the foot warns instead.
    render(
      <ScopePlate
        {...base}
        preview={[]}
        buildDisabled
        footState={{ kind: "warn", text: "Keep at least one value to build" }}
      />,
    );
    expect(screen.queryAllByTestId("ps-seg")).toHaveLength(0);
    expect(
      screen.queryByRole("heading", { level: 2, name: /preview · phase across the series/i }),
    ).not.toBeInTheDocument();
    expect(screen.queryByText(/no clear phase/i)).not.toBeInTheDocument();
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
