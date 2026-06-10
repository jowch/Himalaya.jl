import { forwardRef } from "react";
import type { ButtonHTMLAttributes, ReactNode } from "react";

export type IconButtonTone = "ghost" | "accent" | "danger";

export interface IconButtonProps
  extends Omit<ButtonHTMLAttributes<HTMLButtonElement>, "aria-label"> {
  /** REQUIRED — icon-only buttons have no text; this is the accessible name. */
  label: string;
  /** hover intent. ghost → hover:text-ink (default); accent → hover:text-print-accent
   *  (a rationed terracotta inline action where the brand mark is the point —
   *  NOT close/remove, which is `ghost`/neutral); danger → error-red at rest
   *  (the destructive baseline signal — DESIGN.md §2), strengthening to a subtle
   *  red wash on hover. */
  tone?: IconButtonTone;
  /** render the canonical dismiss glyph (×, U+00D7) as the content. */
  dismiss?: boolean;
  /** glyph or SVG when not a dismiss button (chevron, trash, reorder arrow). */
  children?: ReactNode;
}

const toneClass: Record<IconButtonTone, string> = {
  ghost: "text-ink-faint hover:text-ink",
  accent: "text-ink-faint hover:text-print-accent",
  danger: "text-error hover:bg-error/10",
};

function cx(...parts: Array<string | false | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

// forwardRef so owners can hold the DOM button — e.g. ExportButton returns
// focus to its menu trigger on keyboard close (APG menu-button).
export const IconButton = forwardRef<HTMLButtonElement, IconButtonProps>(function IconButton(
  {
    label,
    tone = "ghost",
    dismiss = false,
    className = "",
    children,
    type = "button",
    title,
    ...props
  },
  ref,
): JSX.Element {
  return (
    <button
      ref={ref}
      type={type}
      aria-label={label}
      title={title ?? label}
      data-tone={tone}
      className={cx(
        // `.icon-button` (styles.css) supplies the >=44px hit-area pseudo-element
        // so the visible box stays compact and never balloons dense chip rows.
        "icon-button relative inline-flex items-center justify-center",
        "rounded p-1 leading-none transition-colors",
        "focus-visible:outline focus-visible:outline-2",
        "focus-visible:outline-offset-2 focus-visible:outline-accent",
        "disabled:cursor-not-allowed disabled:opacity-30",
        toneClass[tone],
        className,
      )}
      {...props}
    >
      {dismiss ? "×" : children}
    </button>
  );
});
