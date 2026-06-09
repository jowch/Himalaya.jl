/**
 * ComposeBar — floating bar that surfaces when ≥1 sample is checked on the
 * contact sheet. Mirrors the CullBar positioning contract: fixed, centred,
 * bottom-20 (sits above CullBar at bottom-7), z-50. The "samples" grain is
 * deliberately different from the CullBar's "frames" grain — never say
 * "frames" here.
 *
 * Presentational contract: NO local state. `count`, `show`, and all handlers
 * are props. The parent page owns selection state.
 *
 * Honest: the bar only renders usefully (aria-hidden absent, buttons tabbable)
 * when `show` is true. `show` MUST equal `count > 0` at the call site.
 *
 * Radius note: uses `rounded-md` (5px) — the single design-system radius step.
 */

import { Button } from "./Button";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface ComposeBarProps {
  /** Number of checked samples. */
  count: number;
  /** true → visible + interactive; false → hidden (faded, pointer-events:none)
   *  while staying mounted so it can animate. Default false. */
  show?: boolean;
  onNewSeries?: () => void;
  onClear?: () => void;
  /** PLACEMENT-ONLY, appended last. */
  className?: string;
}

export function ComposeBar({
  count,
  show = false,
  onNewSeries,
  onClear,
  className,
}: ComposeBarProps): JSX.Element {
  return (
    <div
      data-testid="compose-bar"
      data-show={show ? "true" : "false"}
      aria-hidden={!show || undefined}
      className={cx(
        "fixed left-1/2 bottom-20 -translate-x-1/2 z-50 flex items-center gap-1 bg-ink text-paper rounded-md shadow-lg pl-4 pr-2 py-[7px] transition-opacity",
        show ? "opacity-100" : "opacity-0 pointer-events-none",
        className,
      )}
    >
      <span className="text-sm font-semibold mr-2.5">
        <b className="font-mono">{count}</b>{" "}
        {count === 1 ? "sample" : "samples"}
      </span>
      <Button
        variant="accent"
        onClick={onNewSeries}
        {...(show ? {} : { tabIndex: -1 })}
      >
        + New series
      </Button>
      <Button
        variant="ghostInverse"
        onClick={onClear}
        {...(show ? {} : { tabIndex: -1 })}
      >
        Clear
      </Button>
    </div>
  );
}
