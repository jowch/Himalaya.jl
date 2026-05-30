import { SegmentedControl, type SegmentOption } from "./ui/SegmentedControl";

export type Representation = "waterfall" | "heatmap";

interface RepresentationToggleProps {
  value: Representation;
  onChange: (next: Representation) => void;
}

const OPTIONS: ReadonlyArray<SegmentOption<Representation>> = [
  { value: "waterfall", label: "Waterfall", testId: "repr-waterfall" },
  { value: "heatmap", label: "Heatmap", testId: "repr-heatmap" },
];

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
 *
 * Thin wrapper over the shared SegmentedControl primitive; keeps its
 * `Representation` export + `value`/`onChange` contract so importers are
 * untouched. Active segment is the canonical ink-on-paper fill.
 */
export function RepresentationToggle({ value, onChange }: RepresentationToggleProps): JSX.Element {
  return (
    <SegmentedControl<Representation>
      aria-label="Plot representation"
      role="group"
      variant="bordered"
      testId="representation-toggle"
      options={OPTIONS}
      value={value}
      onChange={onChange}
    />
  );
}
