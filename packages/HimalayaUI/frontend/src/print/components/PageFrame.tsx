import type { ReactNode } from "react";
import { cx } from "../../lib/cx";

export const PAGE_WIDTHS = {
  home: "max-w-[1080px]",
  experiment: "max-w-[1280px]",
} as const;

export type PageWidth = keyof typeof PAGE_WIDTHS;


/** Centers + caps a page body at its surface width. Greenfield pages wrap their
 *  body in this instead of a hand-rolled `mx-auto max-w-[…]`. PLACEMENT-ONLY
 *  className (padding/margin). */
export function PageFrame({
  width,
  className,
  children,
}: {
  width: PageWidth;
  className?: string;
  children: ReactNode;
}): JSX.Element {
  return (
    <div data-testid="page-frame" className={cx("mx-auto w-full", PAGE_WIDTHS[width], className)}>
      {children}
    </div>
  );
}
