import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { Gallery } from "../../src/print/components/Gallery";

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Simple stub child — a plain div with a testid. */
function StubCard({ id }: { id: number }) {
  return <div data-testid={`stub-card-${id}`} />;
}

// ---------------------------------------------------------------------------
// Render with children
// ---------------------------------------------------------------------------

describe("<Gallery> with children", () => {
  it("renders the gallery wall (data-testid=gallery)", () => {
    render(
      <Gallery>
        <StubCard id={1} />
        <StubCard id={2} />
        <StubCard id={3} />
      </Gallery>,
    );
    expect(screen.getByTestId("gallery")).toBeInTheDocument();
  });

  it("renders all 3 passed stub-cards inside the gallery", () => {
    render(
      <Gallery>
        <StubCard id={1} />
        <StubCard id={2} />
        <StubCard id={3} />
      </Gallery>,
    );
    expect(screen.getByTestId("stub-card-1")).toBeInTheDocument();
    expect(screen.getByTestId("stub-card-2")).toBeInTheDocument();
    expect(screen.getByTestId("stub-card-3")).toBeInTheDocument();
  });

  it("does NOT render gallery-empty when children are present", () => {
    render(
      <Gallery>
        <StubCard id={1} />
      </Gallery>,
    );
    expect(screen.queryByTestId("gallery-empty")).not.toBeInTheDocument();
  });
});

// ---------------------------------------------------------------------------
// Empty state: children + empty prop
// ---------------------------------------------------------------------------

describe("<Gallery> with zero children and an empty node", () => {
  it("renders gallery-empty (not the columns wall)", () => {
    render(
      <Gallery empty={<div data-testid="my-empty" />}>
        {/* no children */}
      </Gallery>,
    );
    expect(screen.getByTestId("gallery-empty")).toBeInTheDocument();
  });

  it("renders the provided empty node inside gallery-empty", () => {
    render(
      <Gallery empty={<div data-testid="my-empty" />}>
        {/* no children */}
      </Gallery>,
    );
    expect(screen.getByTestId("my-empty")).toBeInTheDocument();
  });

  it("does NOT render the columns wall (data-testid=gallery) when empty", () => {
    render(
      <Gallery empty={<div data-testid="my-empty" />}>
        {/* no children */}
      </Gallery>,
    );
    expect(screen.queryByTestId("gallery")).not.toBeInTheDocument();
  });

  it("no stub-cards appear in the empty state", () => {
    render(
      <Gallery empty={<div data-testid="my-empty" />}>
        {/* no children */}
      </Gallery>,
    );
    expect(screen.queryByTestId("stub-card-1")).not.toBeInTheDocument();
  });
});

// ---------------------------------------------------------------------------
// Edge case: empty but no `empty` prop provided
// ---------------------------------------------------------------------------

describe("<Gallery> with zero children and no empty prop", () => {
  it("renders the gallery wall (not a crash)", () => {
    render(<Gallery>{/* no children */}</Gallery>);
    expect(screen.getByTestId("gallery")).toBeInTheDocument();
  });

  it("does NOT render gallery-empty", () => {
    render(<Gallery>{/* no children */}</Gallery>);
    expect(screen.queryByTestId("gallery-empty")).not.toBeInTheDocument();
  });
});
