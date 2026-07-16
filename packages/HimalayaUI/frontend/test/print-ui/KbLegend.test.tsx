import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { KbLegend } from "../../src/print/ui/KbLegend";
import type { Shortcut } from "../../src/print/ui/KbLegend";

const shortcuts: Shortcut[] = [
  { keyLabel: "←/→", description: "navigate" },
  { keyLabel: "X", description: "reject" },
  { keyLabel: "Esc", description: "close" },
];

describe("KbLegend", () => {
  it("renders one item per shortcut", () => {
    render(<KbLegend shortcuts={shortcuts} />);
    expect(screen.getAllByTestId("kbkey")).toHaveLength(3);
  });

  it("renders each shortcut's key label and description text", () => {
    render(<KbLegend shortcuts={shortcuts} />);
    for (const s of shortcuts) {
      expect(screen.getByText(s.keyLabel)).toBeTruthy();
      expect(screen.getByText(s.description, { exact: false })).toBeTruthy();
    }
  });

  it('has data-testid="kb-legend"', () => {
    render(<KbLegend shortcuts={shortcuts} />);
    expect(screen.getByTestId("kb-legend")).toBeTruthy();
  });
});
