/**
 * Dock — floating offset control kept reachable while reading full-bleed.
 *
 * Presentational contract: NO local state. `offset` and `onOffsetChange` are
 * props; the parent owns the value. Anchored bottom-right.
 *
 * Layout: an uppercase "Offset" label, a bare `Slider` (no visible label row —
 * its accessible name comes from `ariaLabel="Trace offset"`), and a trailing
 * mono value rendered as `1.20×`.
 */

import { Slider } from "../ui";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface DockProps {
  /** Current trace offset multiplier. Rendered as `${offset.toFixed(2)}×`. */
  offset: number;
  min?: number;
  max?: number;
  step?: number;
  onOffsetChange: (v: number) => void;
  /** PLACEMENT-ONLY, appended last. */
  className?: string;
}

export function Dock({
  offset,
  min = 0.4,
  max = 1.4,
  step = 0.05,
  onOffsetChange,
  className,
}: DockProps): JSX.Element {
  return (
    <div
      data-testid="dock"
      className={cx(
        "fixed right-6 bottom-6 z-50 flex items-center gap-3.5",
        "bg-plate border border-hair-strong rounded-xl px-4 py-3 shadow-lg",
        className,
      )}
    >
      <span className="text-label text-ink-faint">Offset</span>
      <Slider
        value={offset}
        min={min}
        max={max}
        step={step}
        onChange={onOffsetChange}
        ariaLabel="Trace offset"
        className="w-32"
      />
      <span className="text-data font-bold">{`${offset.toFixed(2)}×`}</span>
    </div>
  );
}
