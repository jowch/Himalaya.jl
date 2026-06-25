import type { ReactNode } from "react";
import { Kicker } from "../ui";
import { cx } from "../../lib/cx";


export interface RailSectionProps {
  label: ReactNode;
  count?: ReactNode;
  note?: ReactNode;
  children?: ReactNode;
  className?: string;
}

export function RailSection({ label, count, note, children, className }: RailSectionProps): JSX.Element {
  return (
    <div data-testid="rail-section" className={cx("flex flex-col gap-2.5", className)}>
      <div data-testid="rail-section-head" className="flex items-baseline justify-between">
        <Kicker tone="soft">
          {label}
          {count != null && <span className="text-ink-soft"> {count}</span>}
        </Kicker>
      </div>
      {children}
      {note != null && (
        <div data-testid="rail-section-note" className="text-caption text-ink-soft leading-relaxed">
          {note}
        </div>
      )}
    </div>
  );
}
