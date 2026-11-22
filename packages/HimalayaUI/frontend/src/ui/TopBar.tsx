import type { ReactNode } from "react";
import { cx } from "../lib/cx";

interface TopBarProps {
  wordmark?: ReactNode;
  children?: ReactNode;
  rightSlot?: ReactNode;
  className?: string;
  "data-testid"?: string;
}
export function TopBar({
  wordmark,
  children,
  rightSlot,
  className,
  "data-testid": dataTestid = "topbar",
}: TopBarProps): JSX.Element {
  return (
    <header
      data-testid={dataTestid}
      className={cx(
        "h-14 flex-shrink-0 flex items-center gap-4 px-6 border-b border-hair bg-paper",
        className,
      )}
    >
      {wordmark}
      {children}
      <div className="flex-1" />
      {rightSlot}
    </header>
  );
}