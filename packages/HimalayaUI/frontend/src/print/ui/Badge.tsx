import type { HTMLAttributes, ReactNode } from "react";

interface BadgeProps extends Omit<HTMLAttributes<HTMLSpanElement>, "children"> {
  children: ReactNode;
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** Inline mono count carried alongside a button label (the notes-button "3").
 *  Flat by design: a measured count, not a status pill — no fill, border, or
 *  background. `ink-faint` keeps it quieter than the label it trails. */
export function Badge({
  children,
  className = "",
  ...rest
}: BadgeProps): JSX.Element {
  return (
    <span
      data-testid="badge"
      className={cx("font-mono text-xs text-ink-soft ml-1.5", className)}
      {...rest}
    >
      {children}
    </span>
  );
}
