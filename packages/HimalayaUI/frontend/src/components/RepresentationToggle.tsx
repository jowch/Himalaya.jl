export type Representation = "waterfall" | "heatmap";

interface RepresentationToggleProps {
  value: Representation;
  onChange: (next: Representation) => void;
}

/**
 * Representation segment for the series builder (#175). Waterfall is the only
 * representation the render core (MultiTracePlot) supports today; "Heatmap" is
 * shown disabled because rendering it requires MultiTracePlot/MemberTraceLayer
 * internals, which #175 explicitly carves out to the I3.2 render-core
 * follow-up (#208). The toggle ships so the surface is structurally complete
 * and the heatmap option drops in (enabled) once #208 lands.
 */
export function RepresentationToggle({ value, onChange }: RepresentationToggleProps): JSX.Element {
  return (
    <div
      data-testid="representation-toggle"
      role="radiogroup"
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
        disabled
        title="Heatmap representation is coming soon (#208)"
        className="px-3 py-1.5 text-xs text-ink-faint opacity-50 cursor-not-allowed"
      >
        Heatmap
      </button>
    </div>
  );
}
