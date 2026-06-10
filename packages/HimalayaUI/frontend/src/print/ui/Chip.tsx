import type { ReactNode } from "react";

export type ChipVariant = "static" | "removable" | "add" | "toggle" | "trigger";

/** The chip size scale — orthogonal to {@link ChipVariant}. `sm` is the dense
 *  tag size (text-xs / px-2); `md` is the roomier filter/facet size (text-sm /
 *  px-3). Any variant can be any size. */
export type ChipSize = "sm" | "md";

interface ChipProps {
  variant?: ChipVariant;
  size?: ChipSize;
  children?: ReactNode;
  active?: boolean;
  onClick?: () => void;
  onRemove?: () => void;
  /** Overrides the root element's `data-testid` (default `"chip"`). Lets the
   *  specialized chips that compose this base preserve their own contract
   *  (FilterChip → "filter-chip", FacetChip → "facet-chip", TagPill → "tag-pill"). */
  testId?: string;
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** Shared pill geometry for every variant — the one base that FacetChip,
 *  FilterChip and TagPill collapse onto. Appearance lives here; the consumer's
 *  `className` is placement-only and appended last. Size (padding + text) is an
 *  orthogonal axis applied via {@link sizeClass}, NOT baked into this base. */
const pillBase =
  "inline-flex items-center gap-1 rounded-full font-semibold whitespace-nowrap transition-colors";

/** The size axis — padding + text-size, independent of variant appearance. */
const sizeClass: Record<ChipSize, string> = {
  sm: "px-2 py-px text-xs",
  md: "px-3 py-1 text-sm",
};

/** Base pill in five forms:
 *  - `static`  — a read-only label tag (default).
 *  - `removable` — static + a NEUTRAL trailing × (ink-faint → ink on hover, per
 *    DESIGN.md §2 / review note 9: a remove affordance is not an accent moment).
 *  - `add`     — a dashed "+ add" invite; hover→accent IS correct here (the only
 *    accent on this primitive, the mockup `.tag-add` invitation).
 *  - `toggle`  — a binary selector; active reads as ink/paper INVERSION (not a
 *    hue change), `aria-pressed`/`data-active` reflect state.
 *  - `trigger` — a dropdown opener: `aria-haspopup="menu"` + an aria-hidden ▾.
 *
 *  Variant controls appearance/semantics ONLY; {@link ChipSize} controls
 *  padding + text size. The two axes are independent. */
export function Chip({
  variant = "static",
  size = "md",
  children,
  active = false,
  onClick,
  onRemove,
  testId = "chip",
  className = "",
}: ChipProps): JSX.Element {
  if (variant === "removable") {
    return (
      <span
        data-testid={testId}
        data-variant="removable"
        data-size={size}
        className={cx(
          pillBase,
          sizeClass[size],
          "text-ink-soft bg-plate border border-hair",
          className,
        )}
      >
        {children}
        <button
          type="button"
          aria-label="Remove"
          onClick={onRemove}
          className="text-ink-faint hover:text-ink focus-visible:outline focus-visible:outline-1 focus-visible:outline-offset-0 focus-visible:outline-accent rounded-sm"
        >
          ×
        </button>
      </span>
    );
  }

  if (variant === "add") {
    return (
      <button
        type="button"
        data-testid={testId}
        data-variant="add"
        data-size={size}
        aria-label="Add"
        onClick={onClick}
        className={cx(
          pillBase,
          sizeClass[size],
          "text-ink-faint border border-dashed border-hair-strong hover:text-accent hover:border-accent focus-visible:outline focus-visible:outline-1 focus-visible:outline-offset-0 focus-visible:outline-accent",
          className,
        )}
      >
        {children ?? "+ add"}
      </button>
    );
  }

  if (variant === "toggle") {
    return (
      <button
        type="button"
        data-testid={testId}
        data-variant="toggle"
        data-size={size}
        aria-pressed={active}
        data-active={active ? "true" : "false"}
        onClick={onClick}
        className={cx(
          pillBase,
          sizeClass[size],
          "border focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
          active
            ? "bg-ink text-paper border-ink"
            : "bg-plate text-ink-soft border-hair-strong hover:border-ink-faint",
          className,
        )}
      >
        {children}
      </button>
    );
  }

  if (variant === "trigger") {
    return (
      <button
        type="button"
        data-testid={testId}
        data-variant="trigger"
        data-size={size}
        aria-haspopup="menu"
        onClick={onClick}
        className={cx(
          pillBase,
          sizeClass[size],
          "text-ink bg-plate border border-hair-strong hover:bg-paper-sunk focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent",
          className,
        )}
      >
        {children}
        <span aria-hidden="true" className="text-ink-faint">
          ▾
        </span>
      </button>
    );
  }

  return (
    <span
      data-testid={testId}
      data-variant="static"
      data-size={size}
      className={cx(
        pillBase,
        sizeClass[size],
        "text-ink-soft bg-plate border border-hair",
        className,
      )}
    >
      {children}
    </span>
  );
}
