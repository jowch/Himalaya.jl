import { useId } from "react";
import type { ReactNode } from "react";
import { cx } from "../../lib/cx";

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
  /** Accessible label for the remove button (default `"Remove"`). Thread a
   *  per-tag value such as `"Remove dose 10"` so duplicate pills are
   *  individually named for assistive technology — see LO-TAGDUP. */
  removeLabel?: string;
  /** Overrides the root element's `data-testid` (default `"chip"`). Lets the
   *  specialized chips that compose this base preserve their own contract
   *  (FilterChip → "filter-chip", TagPill → "tag-pill"). */
  testId?: string;
  /** Optional native tooltip text (hover/focus). Used by the `toggle` variant
   *  to explain a filter's meaning; omitted ⇒ no `title` attribute (existing
   *  consumers stay byte-identical). */
  title?: string;
  /** Optional explanatory description (the `toggle` variant only). Rendered as a
   *  visually-hidden span wired via `aria-describedby` so keyboard/SR users hear
   *  the explanation on focus while the accessible NAME stays the label.
   *  Complements `title` (the mouse-only affordance). Omitted ⇒ no
   *  `aria-describedby`/sr-only span (existing consumers stay byte-identical). */
  description?: string;
  className?: string;
}


/** Shared pill geometry for every variant — the one base that FilterChip
 *  and TagPill collapse onto. Appearance lives here; the consumer's
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
  removeLabel = "Remove",
  testId = "chip",
  title,
  description,
  className = "",
}: ChipProps): JSX.Element {
  // Stable id for the optional sr-only description (toggle variant). useId is
  // SSR-safe and collision-free across multiple chips on the page.
  const descId = useId();
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
        {/* WCAG 2.5.8: the × needs a real ≥24×24px hit box, but the pill is a
            dense ~20px tag that must not inflate. h-6/w-6 gives the border box;
            the negative margins collapse its LAYOUT footprint back to roughly
            the bare glyph's (≈8×16px), so pill size is unchanged. -mx-2 is
            tuned to the tightest case (sm, px-2): the button's right edge lands
            exactly at the pill's inner border, never past it — the neighboring
            pill's × and the "+ tag" invite sit only 6px away (gap-1.5), so
            spilling right would overlap another target. Spilling LEFT over the
            pill's own non-interactive text is allowed and intended. */}
        <button
          type="button"
          aria-label={removeLabel}
          onClick={onRemove}
          className="inline-flex h-6 w-6 shrink-0 items-center justify-center -mx-2 -my-1 text-ink-faint hover:text-ink focus-visible:outline focus-visible:outline-1 focus-visible:outline-offset-0 focus-visible:outline-accent rounded-sm"
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
          "text-ink-soft border border-dashed border-hair-strong hover:text-accent hover:border-accent focus-visible:outline focus-visible:outline-1 focus-visible:outline-offset-0 focus-visible:outline-accent",
          className,
        )}
      >
        {children ?? "+ add"}
      </button>
    );
  }

  if (variant === "toggle") {
    return (
      <>
        <button
          type="button"
          data-testid={testId}
          data-variant="toggle"
          data-size={size}
          aria-pressed={active}
          data-active={active ? "true" : "false"}
          {...(title !== undefined ? { title } : {})}
          {...(description !== undefined ? { "aria-describedby": descId } : {})}
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
        {/* Visually-hidden description: keyboard/SR users hear it on focus
            (aria-describedby), the accessible NAME stays the label. `sr-only`
            is the house clip-rect helper — layout-neutral, not appearance. */}
        {description !== undefined && (
          <span id={descId} className="sr-only">
            {description}
          </span>
        )}
      </>
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
