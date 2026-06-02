import type { ReactNode } from "react";
import { Kicker } from "../ui";

export interface PlateHeaderProps {
  /** Optional terracotta eyebrow, e.g. "Integration". */
  kicker?: ReactNode;
  /** Serif title (the plate headline). */
  title: ReactNode;
  /** Optional mono subtitle line. */
  subtitle?: ReactNode;
  /** Heading level for the title (a11y). Default "h2". */
  as?: "h1" | "h2" | "h3";
  /** Right-side tools slot (e.g. a ToolBar). */
  children?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function PlateHeader({
  kicker,
  title,
  subtitle,
  as = "h2",
  children,
  className,
}: PlateHeaderProps): JSX.Element {
  const TitleTag = as;
  return (
    <div
      className={`flex items-start justify-between gap-5 mb-0.5${className ? ` ${className}` : ""}`}
    >
      <div>
        {kicker && (
          <Kicker tone="accent">{kicker}</Kicker>
        )}
        <TitleTag className="text-display">{title}</TitleTag>
        {subtitle && (
          <div className="text-data text-ink-faint mt-1">{subtitle}</div>
        )}
      </div>
      {children && <div>{children}</div>}
    </div>
  );
}
