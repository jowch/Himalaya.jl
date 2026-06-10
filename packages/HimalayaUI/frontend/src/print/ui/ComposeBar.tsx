/**
 * ComposeBar — floating bar that surfaces when ≥1 sample is checked on the
 * contact sheet. Mirrors the CullBar positioning contract: fixed, centred,
 * bottom-20 (sits above CullBar at bottom-7), z-50. The "samples" grain is
 * deliberately different from the CullBar's "frames" grain — never say
 * "frames" here.
 *
 * Presentational contract: NO local state. `count`, `show`, and all handlers
 * are props. The parent page owns selection state. (The ref + layout effect
 * below is the hide mechanism, not state.)
 *
 * Honest: the bar is only interactive (not inert) when `show` is true.
 * `show` MUST equal `count > 0` at the call site.
 *
 * Hidden state: the bar stays mounted to animate, so when `show` is false it
 * is made `inert` (toggled via ref in a layout effect — React 18 has no
 * first-class inert prop). Native inert blurs any focused descendant and
 * removes the subtree from BOTH the a11y tree and the tab order, so a click
 * on Clear that hides the bar cannot strand focus behind an aria-hidden
 * ancestor (the old mechanism's WCAG 4.1.2 failure — see CullBar).
 *
 * Radius note: uses `rounded-md` (5px) — the single design-system radius step.
 */

import { useLayoutEffect, useRef } from "react";
import { Button } from "./Button";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface ComposeBarProps {
  /** Number of checked samples. */
  count: number;
  /** true → visible + interactive; false → hidden (faded, inert,
   *  pointer-events:none) while staying mounted so it can animate.
   *  Default false. */
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
  const barRef = useRef<HTMLDivElement>(null);
  // React 18.3 lacks an inert prop; toggle the attribute before paint instead.
  useLayoutEffect(() => {
    barRef.current?.toggleAttribute("inert", !show);
  }, [show]);

  return (
    <div
      ref={barRef}
      data-testid="compose-bar"
      data-show={show ? "true" : "false"}
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
      <Button variant="accent" onClick={onNewSeries}>
        + New series
      </Button>
      <Button variant="ghostInverse" onClick={onClear}>
        Clear
      </Button>
    </div>
  );
}
