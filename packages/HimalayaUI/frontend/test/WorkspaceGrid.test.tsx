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
  const { container } = render(
    <WorkspaceGrid
      left={<div>L</div>}
      center={<div>C</div>}
      right={<div>R</div>}
    />,
  );
  const center = container.querySelector('[data-slot="center"]')!;
  const right  = container.querySelector('[data-slot="right"]')!;
  const left   = container.querySelector('[data-slot="left"]')!;
  expect(center.className).toMatch(/\border-1\b/);
  expect(right.className).toMatch(/\border-2\b/);
  expect(left.className).toMatch(/\border-3\b/);
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
  expect(container.querySelector('[data-slot="left"]')!.className)
    .toMatch(/\border-1\b/);
  expect(container.querySelector('[data-slot="center"]')!.className)
    .toMatch(/\border-2\b/);
  expect(container.querySelector('[data-slot="right"]')!.className)
    .toMatch(/\border-3\b/);
});

test("uses the new 1400px breakpoint with minmax-floored side columns", () => {
  const { container } = render(
    <WorkspaceGrid left={<div />} center={<div />} right={<div />} />,
  );
  const grid = container.querySelector('[data-testid="workspace-grid"]')!;
  // Side columns floored at 320px so chat / indices stay usable
  expect(grid.className).toContain("min-[1400px]:grid-cols-[minmax(320px,22fr)_56fr_minmax(320px,22fr)]");
});
