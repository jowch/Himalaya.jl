import type { ReactNode } from "react";

/** Single source of truth for per-surface content widths (mockup values).
 *  Changing a surface's width is a one-line edit here. Centred surfaces only —
 *  the focus workspace is full-bleed (work 1fr · rail 350px) and does NOT use a
 *  PageFrame, so it has no entry here. The app bar is full-bleed too. */
export const PAGE_WIDTHS = {
  loupe: "max-w-[1100px]",
  sheet: "max-w-[1240px]",
  folio: "max-w-[1380px]",
  scoping: "max-w-[760px]",
  builder: "max-w-[1180px]",
  home: "max-w-[1080px]",
  experiment: "max-w-[1280px]",
} as const;

export type PageWidth = keyof typeof PAGE_WIDTHS;

function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

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
