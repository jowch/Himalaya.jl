/**
 * CullBar — floating batch-reject action bar.
 *
 * Presentational contract: NO local state. `count`, `show`, and all handlers
 * are props. The parent page owns selection state.
 *
 * Radius note: uses `rounded-md` (5px) — the single design-system radius step —
 * rather than the mockup's `border-radius: 10px`. The collapsed-radius decision
 * is intentional (see docs/design-system.html and the Toast primitive for the
 * same deviation).
 */

import { Button } from "../ui";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface CullBarProps {
  /** Number of currently-selected frames. Shown as the mono count. */
  count: number;
  /** true → visible + interactive; false → hidden (faded, pointer-events:none)
   *  while staying mounted so it can animate. Default false. */
  show?: boolean;
  onReject?: () => void;   // "Drop" (accent)
  onRestore?: () => void;  // "Restore" (ghostInverse)
  onClear?: () => void;    // "Clear" (ghostInverse)
  /** PLACEMENT-ONLY, appended last. */
  className?: string;
}

export function CullBar({
  count,
  show = false,
  onReject,
  onRestore,
  onClear,
  className,
}: CullBarProps): JSX.Element {
  return (
    <div
      data-testid="cull-bar"
      data-show={show ? "true" : "false"}
      aria-hidden={!show}
      className={cx(
        "fixed left-1/2 bottom-7 -translate-x-1/2 z-50 flex items-center gap-1 bg-ink text-paper rounded-md shadow-lg pl-4 pr-2 py-[7px] transition-opacity",
        show ? "opacity-100" : "opacity-0 pointer-events-none",
        className,
      )}
    >
      <span className="text-sm font-semibold mr-2.5">
        <b className="font-mono">{count}</b> frames selected
      </span>
      <Button variant="accent" onClick={onReject} {...(show ? {} : { tabIndex: -1 })}>
        Drop<span className="font-mono text-xs opacity-60 ml-1">X</span>
      </Button>
      <Button variant="ghostInverse" onClick={onRestore} {...(show ? {} : { tabIndex: -1 })}>
        Restore
      </Button>
      <Button variant="ghostInverse" onClick={onClear} {...(show ? {} : { tabIndex: -1 })}>
        Clear<span className="font-mono text-xs opacity-60 ml-1">Esc</span>
      </Button>
    </div>
  );
}
