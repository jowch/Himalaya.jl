import type { ReactNode } from "react";

export function Kbd({ children }: { children: ReactNode }): JSX.Element {
  return (
    <kbd className="inline-flex min-w-[1.5em] whitespace-nowrap justify-center rounded-sm border border-hair-strong border-b-2 bg-plate px-1 text-xs leading-[18px] text-ink-soft">
      {children}
    </kbd>
  );
}
