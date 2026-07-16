import { cx } from "../../lib/cx";
import { KbKey } from "./KbKey";

export interface DockUpLinkProps {
  /** Destination name; rendered as "‹ {label}". */
  label: string;
  onClick: () => void;
  /** When set, render an `<a href>` (with preventDefault) instead of a `<button>`. */
  href?: string;
  /** Key-cap hint (e.g. "esc") shown after the label, so the up-link advertises
   *  its shortcut like every other dock control. */
  kbd?: string;
  /** PLACEMENT-ONLY (e.g. a trailing `mr-1`). Appearance lives here. */
  className?: string;
}

/**
 * The "‹ back" up-link that opens every bottom Dock (Dock grammar §3.3). One
 * appearance, one `data-testid="dock-up-link"` e2e contract; callers vary only
 * the label, the action, the optional key hint, and whether it is a link
 * (`href`) or a button.
 */
export function DockUpLink({ label, onClick, href, kbd, className }: DockUpLinkProps): JSX.Element {
  const cls = cx("inline-flex items-center text-meta font-semibold text-print-accent hover:underline", className);
  const inner = (
    <>
      ‹ {label}
      {kbd && <KbKey className="ml-1.5">{kbd}</KbKey>}
    </>
  );
  if (href !== undefined) {
    return (
      <a
        href={href}
        onClick={(e) => {
          e.preventDefault();
          onClick();
        }}
        className={cls}
        data-testid="dock-up-link"
      >
        {inner}
      </a>
    );
  }
  return (
    <button onClick={onClick} className={cls} data-testid="dock-up-link">
      {inner}
    </button>
  );
}
