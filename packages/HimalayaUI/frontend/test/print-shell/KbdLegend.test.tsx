import { describe, it, expect } from "vitest";
import { render, screen, within } from "@testing-library/react";
import { KbdLegend } from "../../src/print/shell/KbdLegend";
import { SHORTCUTS, shortcutLabel } from "../../src/print/shell/shortcuts";

describe("KbdLegend (registry-driven)", () => {
  it("renders each id's key cap and description straight from the registry", () => {
    render(<KbdLegend ids={["drop", "keep"]} />);
    const legend = screen.getByTestId("kb-legend");
    // descriptions come from SHORTCUTS[id].label
    expect(within(legend).getByText(SHORTCUTS.drop.label)).toBeInTheDocument();
    expect(within(legend).getByText(SHORTCUTS.keep.label)).toBeInTheDocument();
    // key caps come from shortcutLabel(id) — the same string the binding uses
    const caps = within(legend).getAllByTestId("kbkey").map((k) => k.textContent);
    expect(caps).toContain(shortcutLabel("drop"));
    expect(caps).toContain(shortcutLabel("keep"));
  });

  it("forwards placement-only className to the legend row", () => {
    render(<KbdLegend ids={["drop"]} className="mt-4" />);
    expect(screen.getByTestId("kb-legend").className).toContain("mt-4");
  });
});
