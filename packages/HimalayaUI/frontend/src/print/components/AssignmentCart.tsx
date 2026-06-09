import { Children, type ReactNode } from "react";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface AssignmentCartProps {
  /** PhaseBlock children — one per assigned phase. */
  children?: ReactNode;
  /** Override the empty-state copy. */
  empty?: ReactNode;
  /** When provided, shows the "custom index…" footer. */
  onCustomIndex?: () => void;
  /** PLACEMENT-ONLY. Appended last. */
  className?: string;
}

const DEFAULT_EMPTY =
  "No phase assigned — every peak is unindexed. Check a candidate below.";

export function AssignmentCart({
  children,
  empty,
  onCustomIndex,
  className,
}: AssignmentCartProps): JSX.Element {
  const count = Children.count(children);
  return (
    <div
      data-testid="assignment-cart"
      data-phase-count={count}
      className={cx(
        "bg-plate border border-hair rounded overflow-hidden",
        className,
      )}
    >
      {count === 0 ? (
        <div
          data-testid="assignment-empty"
          className="text-caption text-ink-soft p-4"
        >
          {empty ?? DEFAULT_EMPTY}
        </div>
      ) : (
        <>
          {count > 1 && (
            <div
              data-testid="coexistence-tag"
              className="text-kicker text-kicker-faint px-4 pt-2.5"
            >
              Coexistence · {count} phases
            </div>
          )}
          {Children.map(children, (child, i) => (
            <div
              data-testid="cart-block"
              data-divider={i > 0 ? "true" : "false"}
              className={cx("px-4", i > 0 && "border-t border-hair")}
            >
              {child}
            </div>
          ))}
        </>
      )}
      {onCustomIndex && (
        <button
          type="button"
          data-testid="custom-index-trigger"
          onClick={onCustomIndex}
          className="w-full text-left border-t border-dashed border-hair-strong px-4 py-2.5 text-caption font-semibold text-ink-soft transition-colors hover:bg-paper-sunk hover:text-ink"
        >
          <span className="text-print-accent font-bold">+</span> custom index…
        </button>
      )}
    </div>
  );
}
