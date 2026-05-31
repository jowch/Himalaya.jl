import type { ButtonHTMLAttributes, ReactNode } from "react";

export type IconButtonTone = "ghost" | "accent" | "danger";

export interface IconButtonProps
  extends Omit<ButtonHTMLAttributes<HTMLButtonElement>, "aria-label"> {
  /** REQUIRED — icon-only buttons have no text; this is the accessible name. */
  label: string;
  /** hover intent. ghost → hover:text-ink (default); accent → hover:text-print-accent
   *  (a rationed terracotta inline action where the brand mark is the point —
   *  NOT close/remove, which is `ghost`/neutral); danger → hover:text-error (destructive). */
  tone?: IconButtonTone;
  /** render the canonical dismiss glyph (×, U+00D7) as the content. */
  dismiss?: boolean;
  /** glyph or SVG when not a dismiss button (chevron, trash, reorder arrow). */
  children?: ReactNode;
}

const toneClass: Record<IconButtonTone, string> = {
  ghost: "text-ink-faint hover:text-ink",
  accent: "text-ink-faint hover:text-print-accent",
  danger: "text-ink-faint hover:text-error",
};

function cx(...parts: Array<string | false | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export function IconButton({
  label,
  tone = "ghost",
  dismiss = false,
  className = "",
  children,
  type = "button",
  title,
  ...props
}: IconButtonProps): JSX.Element {
  return (
    <button
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
}
