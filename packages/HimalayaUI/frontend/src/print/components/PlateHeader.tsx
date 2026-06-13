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
  /** When the visible `title` is an interactive control (e.g. an editable
   *  Input), pass its plain-text value here. PlateHeader then renders a
   *  visually-hidden heading carrying this text as the real page heading, and
   *  puts the visible `title` in a NON-heading wrapper — so assistive tech gets
   *  a proper, named heading instead of an empty one wrapping a form control.
   *  Omit it (the default) when `title` is static text rendered AS the heading. */
  headingText?: string;
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
  headingText,
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
        {headingText !== undefined ? (
          <>
            {/* The real heading is the title TEXT (sr-only), so heading
                navigation works and no interactive control sits inside a
                heading; the visible interactive `title` rides a plain wrapper. */}
            <TitleTag className="sr-only">{headingText}</TitleTag>
            <div className="text-display">{title}</div>
          </>
        ) : (
          <TitleTag className="text-display">{title}</TitleTag>
        )}
        {subtitle && (
          <div data-role="plate-subtitle" className="text-data text-ink-soft mt-1">
            {subtitle}
          </div>
        )}
      </div>
      {children && <div>{children}</div>}
    </div>
  );
}
