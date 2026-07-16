import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { CombScaffold, MIN_PLOT_W, type ScaffoldRow } from "../../src/print/comb/CombScaffold";

const rows: ScaffoldRow[] = [
  { gutterTitle: "Pn3m", gutterSub: "a = 197 Å" },
  { gutterTitle: "Im3m", gutterSub: "a = 142 Å", preview: true },
];

describe("CombScaffold", () => {
  it("renders one pinned gutter block per row, with the title and sub-line", () => {
    const { container, getByText } = render(
      <CombScaffold rows={rows} xDomain={[0.05, 0.25]} ariaLabel="comb">{() => null}</CombScaffold>,
    );
    expect(container.querySelectorAll('[data-role="gutter-row"]').length).toBe(2);
    expect(getByText("Pn3m")).toBeTruthy();
    expect(getByText("a = 197 Å")).toBeTruthy();
  });

  it("draws a baseline per row; the preview row's baseline is dashed", () => {
    const { container } = render(
      <CombScaffold rows={rows} xDomain={[0.05, 0.25]} ariaLabel="comb">{() => null}</CombScaffold>,
    );
    const baselines = container.querySelectorAll('[data-role="row-baseline"]');
    expect(baselines.length).toBe(2);
    const preview = container.querySelector('[data-role="row-baseline"][data-preview="true"]')!;
    expect(preview.getAttribute("stroke-dasharray")).toBeTruthy();
  });

  it("draws labelled q-axis ticks (major/mid), skipping minor labels", () => {
    const { container } = render(
      <CombScaffold rows={rows} xDomain={[0.05, 0.25]} ariaLabel="comb">{() => null}</CombScaffold>,
    );
    expect(container.querySelectorAll('[data-role="q-tick"]').length).toBeGreaterThan(0);
  });

  it("enforces a min plot width so teeth never crowd (narrow container → wide SVG, pane scrolls)", () => {
    const { container } = render(
      <CombScaffold rows={rows} xDomain={[0.05, 0.25]} maxWidth={320} ariaLabel="comb">{() => null}</CombScaffold>,
    );
    const svg = container.querySelector("svg")!;
    expect(Number(svg.getAttribute("data-plot-w"))).toBeGreaterThanOrEqual(MIN_PLOT_W);
  });

  it("invokes the render-prop with a working log-q projection and a baselineY mapper", () => {
    let seen: { x0: number; x1: number; rowCount: number } | null = null;
    render(
      <CombScaffold rows={rows} xDomain={[0.05, 0.25]} ariaLabel="comb">
        {(ctx) => { seen = { x0: ctx.x.to(0.05), x1: ctx.x.to(0.25), rowCount: ctx.rowCount }; return null; }}
      </CombScaffold>,
    );
    expect(seen!.x1).toBeGreaterThan(seen!.x0); // q increases left→right
    expect(seen!.rowCount).toBe(2);
  });
});
