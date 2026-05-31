import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { EmptyState } from "../../src/print/ui/EmptyState";

describe("EmptyState", () => {
  it("renders the title as a heading", () => {
    render(<EmptyState title="No series match" />);
    expect(screen.getByRole("heading", { name: "No series match" })).toBeTruthy();
  });

  it("renders the body when provided", () => {
    render(
      <EmptyState
        title="No series match"
        body="Try clearing a filter, or widen the beamtime range to see more."
      />,
    );
    expect(
      screen.getByText("Try clearing a filter, or widen the beamtime range to see more."),
    ).toBeTruthy();
  });

  it("does not render a body node when body is omitted", () => {
    render(<EmptyState title="No samples yet" />);
    const root = screen.getByTestId("empty-state");
    // Only the heading is rendered; no body container.
    expect(root.children.length).toBe(1);
    expect(screen.queryByText("anything")).toBeNull();
  });

  it('has data-testid="empty-state"', () => {
    render(<EmptyState title="No samples yet" />);
    expect(screen.getByTestId("empty-state")).toBeTruthy();
  });

  it("forwards a placement className", () => {
    render(<EmptyState title="No samples yet" className="mt-8" />);
    expect(screen.getByTestId("empty-state").className).toContain("mt-8");
  });
});
