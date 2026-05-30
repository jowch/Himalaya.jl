import { SegmentedControl, type SegmentOption } from "./ui/SegmentedControl";

export type ScaleMode = "log" | "linear";

interface ScaleToggleProps {
  value: ScaleMode;
  onChange: (next: ScaleMode) => void;
}

const OPTIONS: ReadonlyArray<SegmentOption<ScaleMode>> = [
  { value: "log", label: "log q", testId: "scale-log" },
  { value: "linear", label: "linear q", testId: "scale-linear" },
];

/**
 * ScaleToggle — log/linear q-axis segmented control (R8 / B-F). Drives
 * MultiTracePlot's `xType`. Thin wrapper over the shared SegmentedControl
 * primitive; keeps its `ScaleMode` export + `value`/`onChange` contract so
 * importers are untouched. Active segment is the canonical ink-on-paper fill.
 */
export function ScaleToggle({ value, onChange }: ScaleToggleProps): JSX.Element {
  return (
    <SegmentedControl<ScaleMode>
      aria-label="q-axis scale"
      role="group"
      variant="bordered"
      testId="scale-toggle"
      options={OPTIONS}
      value={value}
      onChange={onChange}
    />
  );
}
