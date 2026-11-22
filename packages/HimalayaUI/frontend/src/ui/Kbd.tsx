import type { ReactNode } from "react";
import { cx } from "../lib/cx";

export function Kbd({ children, onFill = false }: { children: ReactNode; onFill?: boolean }): JSX.Element {
  return (
    <kbd
      className={cx(
        "inline-flex min-w-[1.5em] whitespace-nowrap justify-center rounded-sm border border-b-2 px-1 text-xs leading-[18px]",
        onFill ? "border-current/40 bg-current/15 text-current" : "border-hair-strong bg-plate text-ink-soft",
      )}
    >
      {children}
    </kbd>
  );
}
