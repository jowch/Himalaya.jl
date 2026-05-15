import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { RowActionZone } from "../src/components/RowActionZone";

describe("RowActionZone — Compare UX E-1", () => {
  it("renders ⋯ overflow and ⋮⋮ drag cue", () => {
    render(<RowActionZone onOverflow={() => {}}/>);
    expect(screen.getByTestId("row-action-overflow")).toBeInTheDocument();
    expect(screen.getByTestId("row-action-drag-cue")).toBeInTheDocument();
  });
  it("dispatches overflow click", () => {
    const onOverflow = vi.fn();
    render(<RowActionZone onOverflow={onOverflow}/>);
    fireEvent.click(screen.getByTestId("row-action-overflow"));
    expect(onOverflow).toHaveBeenCalled();
  });
  it("⋮⋮ is signage (no click handler runs)", () => {
    const onOverflow = vi.fn();
    render(<RowActionZone onOverflow={onOverflow}/>);
    fireEvent.click(screen.getByTestId("row-action-drag-cue"));
    expect(onOverflow).not.toHaveBeenCalled();
  });
});
