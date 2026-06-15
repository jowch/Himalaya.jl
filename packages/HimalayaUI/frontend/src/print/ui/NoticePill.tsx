import type { ReactNode } from "react";

export type NoticePillTone = "new" | "draft";

interface NoticePillProps {
  tone: NoticePillTone;
  children: ReactNode;
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** A compact 10px pill shown in the kick-row of a Series folio card.
 *
 * Two tones:
 *   - `"new"`   — accent-tinted "+N new match" (recipe found new matching samples).
 *   - `"draft"` — dashed, faint "Draft" marker.
 *
 * Appearance (tint, dashed border, accent color) lives here; the consumer's
 * `className` is PLACEMENT-ONLY and appended last — open-placement / closed-look
 * contract.
 *
 * This file is design-guard-EXEMPT (`src/print/ui/`) because it necessarily
 * authors literal appearance (accent tint via `color-mix`, dashed border). */
export function NoticePill({
  tone,
  children,
  className = "",
}: NoticePillProps): JSX.Element {
  const baseClasses =
    "inline-flex items-center rounded-full font-bold tracking-wide whitespace-nowrap";

  if (tone === "new") {
    return (
      <span
        data-testid="notice-pill"
        data-tone="new"
        className={cx(baseClasses, "text-print-accent", className)}
        style={{
          fontSize: "10px",
          letterSpacing: "0.03em",
          padding: "2px 8px",
          background: "color-mix(in oklab, var(--color-accent) 14%, transparent)",
        }}
      >
        {children}
      </span>
    );
  }

  // tone === "draft"
  return (
    <span
      data-testid="notice-pill"
      data-tone="draft"
      className={cx(
        baseClasses,
        "text-ink-soft bg-paper-sunk border border-dashed border-hair-strong",
        className,
      )}
      style={{
        fontSize: "10px",
        letterSpacing: "0.03em",
        padding: "2px 8px",
      }}
    >
      {children}
    </span>
  );
}
