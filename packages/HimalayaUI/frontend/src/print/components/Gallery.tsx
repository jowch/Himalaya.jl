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
  return (
    <div
      data-testid="gallery"
      className={cx("columns-1 sm:columns-2 lg:columns-3 gap-5", className)}
    >
      {React.Children.map(children, (child, i) => (
        <div key={i} className="break-inside-avoid mb-5">
          {child}
        </div>
      ))}
    </div>
  );
}
