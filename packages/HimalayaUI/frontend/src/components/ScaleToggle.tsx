export type ScaleMode = "log" | "linear";

interface ScaleToggleProps {
  value: ScaleMode;
  onChange: (next: ScaleMode) => void;
}

/**
 * ScaleToggle — log/linear q-axis segmented control (R8 / B-F). Drives
 * MultiTracePlot's `xType`. Active segment uses the ink-fill the mockup's
 * `.seg button.on` defines (`background:var(--ink);color:var(--paper)`),
 * matching the sibling RepresentationToggle's active state.
 */
export function ScaleToggle({ value, onChange }: ScaleToggleProps): JSX.Element {
  return (
    <div
      data-testid="scale-toggle"
      role="group"
      aria-label="q-axis scale"
      className="inline-flex overflow-hidden rounded border border-hair-strong"
    >
      <button
        type="button"
        data-testid="scale-log"
        aria-pressed={value === "log"}
        onClick={() => onChange("log")}
        className={`flex-1 px-3 py-1.5 text-xs ${value === "log" ? "bg-ink text-paper" : "text-ink-faint"}`}
      >
        log q
      </button>
      <button
        type="button"
        data-testid="scale-linear"
        aria-pressed={value === "linear"}
        onClick={() => onChange("linear")}
        className={`flex-1 px-3 py-1.5 text-xs ${value === "linear" ? "bg-ink text-paper" : "text-ink-faint"}`}
      >
        linear q
      </button>
    </div>
  );
}
