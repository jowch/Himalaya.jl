/**
 * The grease-pencil reject mark (M-10) — two hand-skewed terracotta strokes
 * drawn over a rejected detector frame, echoing a contact sheet marked up by
 * eye. Ported from sample-table.html's `xmarkSVG`. The stroke colour is the
 * one accent (`--color-print-accent`, terracotta), per DESIGN.md §accent.
 *
 * Pointer-events are off so the mark never intercepts clicks on the thumb.
 */
export function RejectXMark({
  className = "",
  testId = "reject-xmark",
}: {
  className?: string;
  testId?: string;
}): JSX.Element {
  return (
    <svg
      data-testid={testId}
      viewBox="0 0 100 100"
      aria-hidden="true"
      className={`pointer-events-none absolute inset-0 h-full w-full ${className}`}
    >
      <line
        x1="16"
        y1="20"
        x2="86"
        y2="82"
        stroke="var(--color-print-accent)"
        strokeWidth="7"
        strokeLinecap="round"
        transform="rotate(-2 50 50)"
      />
      <line
        x1="84"
        y1="18"
        x2="14"
        y2="84"
        stroke="var(--color-print-accent)"
        strokeWidth="7"
        strokeLinecap="round"
        transform="rotate(2 50 50)"
      />
    </svg>
  );
}
