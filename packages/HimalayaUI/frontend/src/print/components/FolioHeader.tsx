import type { ReactNode } from "react";
import { Kicker } from "../ui";

export interface FolioHeaderProps {
  kicker: string;
  title: string;
  subtitle?: ReactNode;
  count: number;
  countLabel: string;
  /** PLACEMENT-ONLY. Appended last to the root. */
  className?: string;
}

// FIDELITY SNAP: the mockup title is 31px; the named scale has no 31px step
// (and an arbitrary 31px utility is guard-banned), so we use `text-display`
// (27px), the nearest serif role. The count is 26px = `text-headline-lg` exactly.
export function FolioHeader({
  kicker,
  title,
  subtitle,
  count,
  countLabel,
  className,
}: FolioHeaderProps): JSX.Element {
  return (
    <div
      data-testid="folio-header"
      className={`flex items-end justify-between gap-8${className ? ` ${className}` : ""}`}
    >
      <div>
        <Kicker tone="accent">{kicker}</Kicker>
        <h1 className="text-display text-ink">{title}</h1>
        {subtitle && <p className="text-body text-ink-soft mt-2 max-w-[60ch]">{subtitle}</p>}
      </div>
      <div className="shrink-0 text-right">
        <div className="text-headline-lg text-ink leading-none">{count}</div>
        <div className="text-caption text-ink-faint uppercase mt-0.5">{countLabel}</div>
      </div>
    </div>
  );
}
