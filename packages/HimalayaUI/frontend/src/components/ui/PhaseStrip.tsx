import { phaseColor } from "../../phases";

/**
 * PhaseStrip — one cell per series member coloured by that member's confirmed
 * phase, captioned with the phase story (single "throughout", "A → B"
 * transition, or an empty label). The shared primitive behind the Series folio
 * card strip and the series-scoping preview strip.
 *
 * Closed look / open placement: appearance (size, colours, gradient) is owned
 * here; the consumer `className` is PLACEMENT ONLY (margin / width / grid).
 * Unindexed cells use the canonical `var(--color-ink-faint)` tone (matches the
 * mini-waterfall's UNINDEXED_COLOR). The "throughout vs transition" caption is
 * derived from the COUNT OF DISTINCT indexed phases (the truthful rule): a
 * non-monotone strip like [Pn3m, Lamellar, Pn3m] reads as a transition.
 *
 * Per-segment `aria-label`/`title` carry the phase name as the accessible
 * second channel (colour is never the sole signal — see phases.ts). The visual
 * glyph/pattern second channel is deferred to the plotting redesign.
 */

/** One per series member, in display (low→high) order. `coexistWith` drives the
 *  2-stop gradient; a null phase = unindexed (the neutral ink-faint cell). */
export interface PhaseSegment {
  phase: string | null;
  coexistWith?: string | null;
  /** Durable 3-state assignment (Plan E E-5/E-7). `form_factor` renders a hollow
   *  dashed cell (a real trace but no Bragg peaks); `null` renders a faint,
   *  distinct cell. Omitted → a plain indexed / unindexed cell. */
  state?: "form_factor" | "null";
}

export type PhaseStripSize = "sm" | "md"; // sm = legacy 8px Scoping bar; md = 7px folio bar (default)

export interface PhaseStripProps {
  segments: PhaseSegment[];
  /** Discrete size. "md" (7px, default) | "sm" (8px legacy Scoping bar). */
  size?: PhaseStripSize;
  /** Caption when no segment is indexed. Default "No clear phase". */
  emptyLabel?: string;
  /** PLACEMENT ONLY: margin / width / grid position. No appearance utilities. */
  className?: string;
}

const UNINDEXED = "var(--color-ink-faint)";

function cx(...parts: Array<string | false | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

const sizeClass: Record<PhaseStripSize, string> = {
  sm: "h-2 gap-0.5", // 8px bar
  md: "h-[7px] gap-[2px]", // 7px bar (folio canonical)
};

function segBackground(seg: PhaseSegment): string {
  if (seg.phase === null) return UNINDEXED;
  if (seg.coexistWith) {
    return `linear-gradient(100deg, ${phaseColor(seg.phase)} 42%, ${phaseColor(seg.coexistWith)} 58%)`;
  }
  return phaseColor(seg.phase);
}

function segLabel(seg: PhaseSegment): string {
  if (seg.state === "form_factor") return "Form factor (no Bragg peaks)";
  if (seg.state === "null") return "No phase";
  if (seg.phase === null) return "Unindexed";
  if (seg.coexistWith) return `${seg.phase} + ${seg.coexistWith} (coexistence)`;
  return seg.phase;
}

export function PhaseStrip({
  segments,
  size = "md",
  emptyLabel = "No clear phase",
  className = "",
}: PhaseStripProps): JSX.Element {
  const indexed = segments
    .map((s) => s.phase)
    .filter((p): p is string => p !== null);
  const first = indexed.length > 0 ? indexed[0]! : null;
  const last = indexed.length > 0 ? indexed[indexed.length - 1]! : null;
  const distinct = new Set(indexed);

  return (
    <div className={className} data-size={size}>
      <div className={cx("flex", sizeClass[size])}>
        {segments.map((seg, i) => {
          // Form-factor → hollow dashed cell; null → a faint distinct cell.
          if (seg.state === "form_factor") {
            return (
              <div
                key={i}
                data-testid="ps-seg"
                data-state="form_factor"
                aria-label={segLabel(seg)}
                title={segLabel(seg)}
                className="flex-1 rounded-[1.5px] border border-dashed border-hair-strong bg-transparent"
              />
            );
          }
          if (seg.state === "null") {
            return (
              <div
                key={i}
                data-testid="ps-seg"
                data-state="null"
                aria-label={segLabel(seg)}
                title={segLabel(seg)}
                className="flex-1 rounded-[1.5px] bg-hair"
              />
            );
          }
          return (
            <div
              key={i}
              data-testid="ps-seg"
              aria-label={segLabel(seg)}
              title={segLabel(seg)}
              className="flex-1 rounded-[1.5px]"
              style={{ background: segBackground(seg) }}
            />
          );
        })}
      </div>
      <div
        data-testid="ps-cap"
        className="mt-1.5 flex items-center gap-1.5 text-base text-ink-soft"
      >
        {first === null || last === null ? (
          <span className="font-semibold text-ink-faint">{emptyLabel}</span>
        ) : distinct.size > 1 ? (
          <>
            <span className="font-semibold" style={{ color: phaseColor(first) }}>
              {first}
            </span>
            <span className="text-ink-faint" aria-hidden="true">
              →
            </span>
            <span className="font-semibold" style={{ color: phaseColor(last) }}>
              {last}
            </span>
          </>
        ) : (
          <>
            <span className="font-semibold" style={{ color: phaseColor(first) }}>
              {first}
            </span>
            <span className="text-ink-faint">throughout</span>
          </>
        )}
      </div>
    </div>
  );
}
