export type Representation = "waterfall" | "heatmap";

interface RepresentationToggleProps {
  value: Representation;
  onChange: (next: Representation) => void;
}

/**
 * Representation segment for the series builder. Picks the plot's layout
 * vocabulary in `MultiTracePlot` (#208):
 *
 *   - `waterfall` — stacked 1-D traces (the legacy, the publication figure)
 *   - `heatmap`   — q-binned intensity rows (peaks-only; surfaces migration)
 *
 * "Coloring vs layout" stays separate: this is the layout axis; coloring
 * lives in `GroupingModeToggle`. They compose: a heatmap can be byPhase or
 * bySample, a waterfall can be either.
 *
 * Originally shipped with the heatmap button disabled (#175 carved out the
 * render-core work to #208); both modes are now live.
 */
export function RepresentationToggle({ value, onChange }: RepresentationToggleProps): JSX.Element {
  return (
    <div
      data-testid="representation-toggle"
      role="group"
      aria-label="Plot representation"
      className="inline-flex overflow-hidden rounded border border-hair-strong"
    >
      <button
        type="button"
        data-testid="repr-waterfall"
        aria-pressed={value === "waterfall"}
        onClick={() => onChange("waterfall")}
        className={`px-3 py-1.5 text-xs ${value === "waterfall" ? "bg-ink text-paper" : "text-ink-faint"}`}
      >
        Waterfall
      </button>
      <button
        type="button"
        data-testid="repr-heatmap"
        aria-pressed={value === "heatmap"}
        onClick={() => onChange("heatmap")}
        className={`px-3 py-1.5 text-xs ${value === "heatmap" ? "bg-ink text-paper" : "text-ink-faint"}`}
      >
        Heatmap
      </button>
    </div>
  );
}
