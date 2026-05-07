/**
 * Observable Plot scale-cast helpers (extracted in Phase 6, Task 6.2).
 *
 * The plot element returned by `Plot.plot(...)` carries a runtime `.scale(name)`
 * method that returns an object with `invert(px)` and `apply(value)` callables.
 * These are NOT in `@types/observablehq__plot` (or the DOM types), so callers
 * must use an `unknown` cast — see CLAUDE.md frontend gotchas
 * ("Observable Plot inside React").
 *
 * `invertQ(plot, px)` returns the `q` value for a pixel x-coordinate, or
 * `null` when the plot is missing or the scale isn't present yet (the plot
 * element can be transiently null between render and ResizeObserver tick).
 *
 * `applyQ(plot, q)` is the forward direction — pixel for a q value — and is
 * a thin wrapper for symmetry. Both are pure and side-effect-free.
 */

type ScaleLike = {
  invert?: (px: number) => number;
  apply?: (value: number) => number;
};

type PlotLike = {
  scale: (name: string) => ScaleLike | undefined;
};

function asPlot(plot: unknown): PlotLike | null {
  if (plot === null || plot === undefined) return null;
  const candidate = plot as { scale?: unknown };
  if (typeof candidate.scale !== "function") return null;
  return plot as PlotLike;
}

/**
 * Pixel → q. Returns null when the plot or its x-scale is unavailable, or
 * when the input pixel isn't finite.
 */
export function invertQ(plot: unknown, px: number): number | null {
  if (!Number.isFinite(px)) return null;
  const p = asPlot(plot);
  if (!p) return null;
  const scale = p.scale("x");
  if (!scale || typeof scale.invert !== "function") return null;
  return scale.invert(px);
}

/**
 * q → pixel. Mirror of `invertQ`. Returns null when the scale isn't ready.
 */
export function applyQ(plot: unknown, q: number): number | null {
  if (!Number.isFinite(q)) return null;
  const p = asPlot(plot);
  if (!p) return null;
  const scale = p.scale("x");
  if (!scale || typeof scale.apply !== "function") return null;
  return scale.apply(q);
}
