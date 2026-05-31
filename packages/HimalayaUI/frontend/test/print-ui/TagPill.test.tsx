import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { TagPill } from "../../src/print/ui/TagPill";

describe("TagPill", () => {
  it("renders its children text", () => {
    render(<TagPill>lipid-A</TagPill>);
    expect(screen.getByTestId("tag-pill").textContent).toContain("lipid-A");
  });

  it('has data-testid="tag-pill"', () => {
    render(<TagPill>lipid-A</TagPill>);
    expect(screen.getByTestId("tag-pill")).toBeTruthy();
  });

  it("renders NO remove button when onRemove is omitted", () => {
    render(<TagPill>lipid-A</TagPill>);
    expect(screen.queryByRole("button", { name: "Remove tag" })).toBeNull();
  });

  it("renders a remove button when onRemove is provided and fires it on click", () => {
    const onRemove = vi.fn();
    render(<TagPill onRemove={onRemove}>lipid-A</TagPill>);
    const button = screen.getByRole("button", { name: "Remove tag" });
    expect(button).toBeTruthy();
    fireEvent.click(button);
    expect(onRemove).toHaveBeenCalledTimes(1);
  });
});
