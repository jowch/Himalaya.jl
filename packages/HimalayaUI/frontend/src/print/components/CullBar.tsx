/**
 * CullBar — floating batch-verdict action bar (Drop / Keep / Restore).
 *
 * Presentational contract: NO local state. `count`, `show`, and all handlers
 * are props. The parent page owns selection state. (The ref + layout effect
 * below is the hide mechanism, not state.)
 *
 * Hidden state: the bar stays mounted to animate, so when `show` is false it
 * is made `inert` (toggled via ref in a layout effect — React 18 has no
 * first-class inert prop). Native inert blurs any focused descendant and
 * removes the subtree from BOTH the a11y tree and the tab order in one move.
 * The previous aria-hidden + per-button tabIndex=-1 pair failed WCAG 4.1.2:
 * clicking Drop/Restore/Clear empties the selection, `show` flips false, and
 * aria-hidden landed on an element whose descendant still held focus (Chrome
 * blocks it and logs a warning). inert handles that ordering itself.
 *
 * Radius note: uses `rounded-md` (5px) — the single design-system radius step —
 * rather than the mockup's `border-radius: 10px`. The collapsed-radius decision
 * is intentional (see docs/design-system.html and the Toast primitive for the
 * same deviation).
 */

import { useLayoutEffect, useRef } from "react";
import { Button } from "../ui";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface CullBarProps {
  /** Number of currently-selected frames. Shown as the mono count. */
  count: number;
  /** Number of distinct samples the selection spans (SA-F3). When > 1 the
   *  count line discloses the spread ("… across M samples") so a selection
   *  accumulated across off-screen rows is never silently included in a
   *  Drop. Omitted or 1 → no suffix: no noise when nothing is hidden. The
   *  parent must derive this from the SAME mapping its Drop handler uses. */
  sampleCount?: number;
  /** true → visible + interactive; false → hidden (faded, inert,
   *  pointer-events:none) while staying mounted so it can animate.
   *  Default false. */
  show?: boolean;
  onReject?: () => void;   // "Drop" (accent)
  onKeep?: () => void;     // "Keep" (success — constructive, never a second accent)
  onRestore?: () => void;  // "Restore" (ghostInverse)
  onClear?: () => void;    // "Clear" (ghostInverse)
  /** PLACEMENT-ONLY, appended last. */
  className?: string;
}

export function CullBar({
  count,
  sampleCount,
  show = false,
  onReject,
  onKeep,
  onRestore,
  onClear,
  className,
}: CullBarProps): JSX.Element {
  const barRef = useRef<HTMLDivElement>(null);
  // React 18.3 lacks an inert prop; toggle the attribute before paint instead.
  useLayoutEffect(() => {
    barRef.current?.toggleAttribute("inert", !show);
  }, [show]);

  return (
    <div
      ref={barRef}
      data-testid="cull-bar"
      data-show={show ? "true" : "false"}
      className={cx(
        "fixed left-1/2 bottom-7 -translate-x-1/2 z-50 flex items-center gap-1 bg-ink text-paper rounded-md shadow-lg pl-4 pr-2 py-[7px] transition-opacity",
        show ? "opacity-100" : "opacity-0 pointer-events-none",
        className,
      )}
    >
      <span className="text-sm font-semibold mr-2.5">
        <b className="font-mono">{count}</b> frame{count === 1 ? "" : "s"} selected
        {sampleCount !== undefined && sampleCount > 1 && (
          // Spread disclosure (SA-F3). Rendered only for a spread (> 1), so
          // "samples" is always correctly plural here — "across 1 samples"
          // can never appear.
          <>
            {" across "}
            <b className="font-mono">{sampleCount}</b>
            {" samples"}
          </>
        )}
      </span>
      <Button variant="accent" onClick={onReject}>
        Drop<span className="font-mono text-xs opacity-60 ml-1">X</span>
      </Button>
      <Button variant="success" onClick={onKeep}>
        Keep<span className="font-mono text-xs opacity-60 ml-1">K</span>
      </Button>
      <Button variant="ghostInverse" onClick={onRestore}>
        Restore
      </Button>
      <Button variant="ghostInverse" onClick={onClear}>
        Clear<span className="font-mono text-xs opacity-60 ml-1">Esc</span>
      </Button>
    </div>
  );
}
