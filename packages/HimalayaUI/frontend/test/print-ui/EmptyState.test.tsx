import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { EmptyState } from "../../src/print/ui/EmptyState";

describe("EmptyState", () => {
  it("renders the title as a heading", () => {
    render(<EmptyState title="No series match" />);
    expect(screen.getByRole("heading", { name: "No series match" })).toBeTruthy();
  });

  it("defaults the title to an h2 (EmptyState usually nests under a page h1)", () => {
    render(<EmptyState title="No series match" />);
    expect(screen.getByRole("heading", { level: 2, name: "No series match" })).toBeTruthy();
  });

  it("renders the title as an h1 when as=\"h1\" (EmptyState IS the page, e.g. not-found)", () => {
    // FO-NFHEAD: on a route where EmptyState is the only heading, an h2 leaves
    // the document with no h1 (a level-skip). The consumer raises it to h1.
    render(<EmptyState title="Sample not found" as="h1" />);
    expect(screen.getByRole("heading", { level: 1, name: "Sample not found" })).toBeTruthy();
    expect(screen.queryByRole("heading", { level: 2 })).toBeNull();
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

  it("renders the action slot below the body", () => {
    render(
      <EmptyState
        title="Could not load the corpus"
        body="The sample list request failed."
        action={<button>Reload the corpus</button>}
      />,
    );
    const root = screen.getByTestId("empty-state");
    const button = screen.getByRole("button", { name: "Reload the corpus" });
    expect(root.contains(button)).toBe(true);
    // The action container is the last child (below the body).
    expect(root.lastElementChild!.contains(button)).toBe(true);
    expect(root.lastElementChild!.textContent).not.toContain(
      "The sample list request failed.",
    );
  });

  it("does not render an action container when action is omitted", () => {
    render(<EmptyState title="No samples yet" body="Nothing here." />);
    const root = screen.getByTestId("empty-state");
    // heading + body only
    expect(root.children.length).toBe(2);
  });
});
