import { KbKey } from "./KbKey";

export interface Shortcut {
  keyLabel: string;
  description: string;
}

interface KbLegendProps {
  shortcuts: Shortcut[];
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** KbLegend — a flat row of shortcut hints, each pairing a KbKey cap with a
 *  short description. Composes KbKey (DRY); the legend is a single flat row of
 *  tonal/voice separation, no plate (checklist H).
 *
 *  The legend itself is non-interactive (static reference of available
 *  shortcuts), so checklist C (interaction states) and D (touch target) are N/A.
 *
 *  Voice (checklist E): key caps are mono via KbKey (measured/literal tokens);
 *  descriptions are sans prose chrome — the two correct, distinct voices. */
export function KbLegend({ shortcuts, className = "" }: KbLegendProps): JSX.Element {
  return (
    <div
      data-testid="kb-legend"
      className={cx("flex flex-wrap gap-5 text-sm text-ink-faint", className)}
    >
      {shortcuts.map((s) => (
        <span key={s.keyLabel} className="inline-flex items-center gap-1.5">
          <KbKey>{s.keyLabel}</KbKey>
          {s.description}
        </span>
      ))}
    </div>
  );
}
