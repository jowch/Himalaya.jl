import { cx } from "../../lib/cx";

export interface DockUpLinkProps {
  /** Destination name; rendered as "‹ {label}". */
  label: string;
  onClick: () => void;
  /** When set, render an `<a href>` (with preventDefault) instead of a `<button>`. */
  href?: string;
  /** PLACEMENT-ONLY (e.g. a trailing `mr-1`). Appearance lives here. */
  className?: string;
}

/**
 * The "‹ back" up-link that opens every bottom Dock (Dock grammar §3.3). One
 * appearance, one `data-testid="dock-up-link"` e2e contract; callers vary only
 * the label, the action, and whether it is a link (`href`) or a button.
 */
export function DockUpLink({ label, onClick, href, className }: DockUpLinkProps): JSX.Element {
  const cls = cx("text-meta font-semibold text-print-accent hover:underline", className);
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
        ‹ {label}
      </a>
    );
  }
  return (
    <button onClick={onClick} className={cls} data-testid="dock-up-link">
      ‹ {label}
    </button>
  );
}
