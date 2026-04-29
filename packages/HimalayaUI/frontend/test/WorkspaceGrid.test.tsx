import { render, screen } from "@testing-library/react";
import { WorkspaceGrid } from "../src/components/WorkspaceGrid";

test("renders the three slots in their slots", () => {
  render(
    <WorkspaceGrid
      left={<div data-testid="L">L</div>}
      center={<div data-testid="C">C</div>}
      right={<div data-testid="R">R</div>}
    />,
  );
  expect(screen.getByTestId("L")).toBeInTheDocument();
  expect(screen.getByTestId("C")).toBeInTheDocument();
  expect(screen.getByTestId("R")).toBeInTheDocument();

  // Slots are tagged with data-slot so consumers can rely on them in E2E if needed
  expect(screen.getByTestId("workspace-grid").querySelector('[data-slot="left"]'))
    .toContainElement(screen.getByTestId("L"));
  expect(screen.getByTestId("workspace-grid").querySelector('[data-slot="center"]'))
    .toContainElement(screen.getByTestId("C"));
  expect(screen.getByTestId("workspace-grid").querySelector('[data-slot="right"]'))
    .toContainElement(screen.getByTestId("R"));
});

test("default mobileOrder puts center first, right second, left last", () => {
  // We assert on the data-mobile-order attribute (not on Tailwind class
  // strings) so the test stays meaningful even if the styling vocabulary
  // changes. The component still emits the corresponding order-N Tailwind
  // class — that's what actually drives the layout — but tests should not
  // couple to its name. See CLAUDE.md: "Never assert on Tailwind class
  // strings — they change when styling evolves."
  const { container } = render(
    <WorkspaceGrid
      left={<div>L</div>}
      center={<div>C</div>}
      right={<div>R</div>}
    />,
  );
  expect(container.querySelector('[data-slot="center"]'))
    .toHaveAttribute("data-mobile-order", "1");
  expect(container.querySelector('[data-slot="right"]'))
    .toHaveAttribute("data-mobile-order", "2");
  expect(container.querySelector('[data-slot="left"]'))
    .toHaveAttribute("data-mobile-order", "3");
});

test("custom mobileOrder reorders the slots", () => {
  const { container } = render(
    <WorkspaceGrid
      left={<div>L</div>}
      center={<div>C</div>}
      right={<div>R</div>}
      mobileOrder={["left", "center", "right"]}
    />,
  );
  expect(container.querySelector('[data-slot="left"]'))
    .toHaveAttribute("data-mobile-order", "1");
  expect(container.querySelector('[data-slot="center"]'))
    .toHaveAttribute("data-mobile-order", "2");
  expect(container.querySelector('[data-slot="right"]'))
    .toHaveAttribute("data-mobile-order", "3");
});

// Note: the previous "uses the new 1400px breakpoint with minmax-floored
// side columns" test was removed. Asserting that a specific
// `min-[1400px]:grid-cols-[...]` Tailwind class string is present in the
// className is testing implementation, not behavior. The breakpoint
// behavior is verified by manual + E2E review at multiple viewport widths.
// If we ever want a programmatic gate, a Playwright test that resizes and
// asserts on layout (computed grid template, slot counts) is the right
// shape — the className grep is brittle either way.
