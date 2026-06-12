import React from "react";
import type { ReactNode } from "react";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface GalleryProps {
  /** SeriesCard elements, already filtered/sorted by the page. Gallery only lays them out. */
  children?: ReactNode;
  /** Rendered (instead of the wall) when there are zero children. */
  empty?: ReactNode;
  /** PLACEMENT-ONLY, appended last. */
  className?: string;
}

/**
 * The Series-folio masonry wall: a CSS multi-column layout of SeriesCard
 * children with distinct heights (a 6-sample card is taller than a 3-sample
 * one, so the wall reads as a masonry of real figures, not a uniform grid).
 * When the (already-filtered) list is empty and an `empty` node is provided,
 * it renders that instead of the columns wall.
 *
 * Column gap: `gap-5` (20px) — the mockup specifies 22px; nearest named step.
 * Child wrapper bottom margin: `mb-5` (20px) — matches the column gap so the
 * gutter is visually uniform in both directions.
 *
 * Design-guard contract: `src/print/components/**` is NOT guard-exempt.
 * Only named-role classes and named layout utilities are used here.
 */
export function Gallery({
  children,
  empty,
  className,
}: GalleryProps): JSX.Element {
  const isEmpty = React.Children.count(children) === 0;

  // Empty state: only shown when there are zero children AND an empty node is
  // provided. The columns wall is not rendered in this branch.
  if (isEmpty && empty != null) {
    return <div data-testid="gallery-empty">{empty}</div>;
  }

  // Normal: lay out children in a responsive masonry-style multi-column grid.
  // Breakpoints: 1 col (default) → 2 col (sm ≥ 640px) → 3 col (lg ≥ 1024px).
  // Each child is wrapped in break-inside-avoid + mb-5 so cards don't split
  // across columns and the bottom gutter matches the column gap.
  //
  // Sparse-wall reflow: cap the effective column count to the number of cards.
  // The full 3-column wall assumes a dense corpus; with 1–2 cards the surplus
  // columns stay empty and the wall reads as broken (everything piled left).
  // Capping `lg:columns-N` (and `sm:columns-N`) to the card count keeps a small
  // corpus balanced — a 2-card folio fills 2 columns, not 1-of-3.
  const count = React.Children.count(children);
  const colClass = columnClassFor(count);

  return (
    <div
      data-testid="gallery"
      data-columns={Math.min(count, 3)}
      className={cx(colClass, "gap-5", className)}
    >
      {React.Children.map(children, (child, i) => (
        <div key={i} className="break-inside-avoid mb-5">
          {child}
        </div>
      ))}
    </div>
  );
}

/** Responsive column track classes, capped to the card count so a sparse wall
 *  stays balanced. 1 card → single column at every width; 2 cards → up to 2;
 *  3+ cards → the full 1→2→3 ladder. Static class strings (no interpolation) so
 *  Tailwind's content scan keeps each utility. */
function columnClassFor(count: number): string {
  if (count <= 1) return "columns-1";
  if (count === 2) return "columns-1 sm:columns-2";
  return "columns-1 sm:columns-2 lg:columns-3";
}
