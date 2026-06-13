import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { useEffect } from "react";
import { Gallery } from "../../src/print/components/Gallery";

describe("<Gallery> reorder identity (FOL-RESCORE5)", () => {
  it("preserves card identity across a re-sort (wrapper key follows the child, not its position)", () => {
    let mounts = 0;
    function Counted({ id }: { id: number }) {
      useEffect(() => {
        mounts += 1;
      }, []);
      return <div data-testid={`counted-${id}`} />;
    }
    const { rerender } = render(
      <Gallery>
        {[<Counted key="a" id={1} />, <Counted key="b" id={2} />]}
      </Gallery>,
    );
    expect(mounts).toBe(2);
    // Re-sort: swap the two cards. With positional wrapper keys the children
    // would remount in place (mounts -> 4); with child-derived keys the whole
    // wrapper moves and the card subtree is preserved (mounts stays 2).
    rerender(
      <Gallery>
        {[<Counted key="b" id={2} />, <Counted key="a" id={1} />]}
      </Gallery>,
    );
    expect(mounts).toBe(2);
  });
});

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
// Sparse-wall column cap (F3): the effective column count is capped to the
// number of cards so a 1–2 card wall stays balanced (data-columns reflects it).
// ---------------------------------------------------------------------------

describe("<Gallery> sparse-wall column cap", () => {
  it("a single card caps the wall to one column (data-columns=1)", () => {
    render(
      <Gallery>
        <StubCard id={1} />
      </Gallery>,
    );
    expect(screen.getByTestId("gallery")).toHaveAttribute("data-columns", "1");
  });

  it("two cards cap the wall to two columns (data-columns=2)", () => {
    render(
      <Gallery>
        <StubCard id={1} />
        <StubCard id={2} />
      </Gallery>,
    );
    expect(screen.getByTestId("gallery")).toHaveAttribute("data-columns", "2");
  });

  it("three or more cards use the full three-column wall (data-columns=3)", () => {
    render(
      <Gallery>
        <StubCard id={1} />
        <StubCard id={2} />
        <StubCard id={3} />
        <StubCard id={4} />
      </Gallery>,
    );
    expect(screen.getByTestId("gallery")).toHaveAttribute("data-columns", "3");
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
