import type { ReactNode } from "react";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

interface TopBarProps {
  wordmark?: ReactNode;
  children?: ReactNode;
  rightSlot?: ReactNode;
  className?: string;
}

/**
 * Fixed 56px (`h-14`) header shell: flat, paper background, bottom hairline,
 * never shrinks. Layout order is wordmark → `children` (stage tabs / facet /
 * stepper) → a `flex-1` spacer that pushes `rightSlot` to the far right. The
 * shell itself is non-interactive — its slot contents own their interaction.
 */
export function TopBar({
  wordmark,
  children,
  rightSlot,
  className,
}: TopBarProps): JSX.Element {
  return (
    <header
      data-testid="topbar"
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
