import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { RowActionZone } from "../src/components/RowActionZone";

describe("RowActionZone", () => {
  it("renders the ⋮⋮ drag cue as inert signage", () => {
    render(<RowActionZone />);
    const cue = screen.getByTestId("row-action-drag-cue");
    expect(cue).toBeInTheDocument();
    expect(cue).toHaveAttribute("aria-hidden", "true");
  });
  // M2: the ⋯ overflow button toggled a `data-overflow-open` attribute that
  // nothing consumed (opened no menu) — a no-op affordance. Removed.
  it("no longer renders the no-op ⋯ overflow button", () => {
    render(<RowActionZone />);
    expect(screen.queryByTestId("row-action-overflow")).toBeNull();
  });
});
