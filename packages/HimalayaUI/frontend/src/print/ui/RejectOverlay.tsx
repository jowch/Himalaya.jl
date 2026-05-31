function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

interface RejectOverlayProps {
  className?: string;
}

/** The hand-skewed terracotta grease-pencil ✕ over a rejected detector frame.
 *  Decorative only: it scales to fill its `relative` parent and is `aria-hidden`
 *  because the rejected STATE is announced elsewhere (e.g. a reason chip / label).
 *  The consumer applies the frame dimming; rejection thus reads via the ✕ SHAPE
 *  + the dimmed frame, never hue alone (checklist A). Accent here is the rationed
 *  reject grease-pencil mark (checklist B). */
export function RejectOverlay({ className = "" }: RejectOverlayProps): JSX.Element {
  return (
    <svg
      data-testid="reject-overlay"
      aria-hidden="true"
      viewBox="0 0 100 100"
      preserveAspectRatio="none"
      className={cx("absolute inset-0 w-full h-full pointer-events-none", className)}
    >
      <line x1="16" y1="20" x2="86" y2="82" stroke="var(--color-accent)" strokeWidth={7} strokeLinecap="round" transform="rotate(-2 50 50)" />
      <line x1="84" y1="18" x2="14" y2="84" stroke="var(--color-accent)" strokeWidth={7} strokeLinecap="round" transform="rotate(2 50 50)" />
    </svg>
  );
}
