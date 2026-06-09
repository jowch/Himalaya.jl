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
 *  N-band coexistence gradient; a null phase = unindexed (the neutral
 *  ink-faint cell). */
export interface PhaseSegment {
  phase: string | null;
  /** Additional coexisting phases beyond the dominant `phase` (in addition
   *  order). Drives an N-stop equal-band gradient; the dominant `phase` is the
   *  first band. Empty/omitted/null = single-phase cell. Supports 2-, 3-,
   *  N-phase coexistence. */
  coexistWith?: string[] | null;
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
  /** Orientation (Plan E E-5). `"horizontal"` (default) lays cells left→right
   *  with the caption beneath (folio bar); `"vertical"` stacks cells top→bottom
   *  (the Series waterfall companion) and drops the caption. */
  orientation?: "horizontal" | "vertical";
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
  const coexist = seg.coexistWith ?? [];
  if (coexist.length > 0) {
    const all = [seg.phase, ...coexist];
    const n = all.length;
    const stops = all
      .map((p, i) => {
        const start = ((i / n) * 100).toFixed(2);
        const end = (((i + 1) / n) * 100).toFixed(2);
        return `${phaseColor(p)} ${start}%, ${phaseColor(p)} ${end}%`;
      })
      .join(", ");
    return `linear-gradient(100deg, ${stops})`;
  }
  return phaseColor(seg.phase);
}

function segLabel(seg: PhaseSegment): string {
  if (seg.state === "form_factor") return "Form factor (no Bragg peaks)";
  if (seg.state === "null") return "No phase";
  if (seg.phase === null) return "Unindexed";
  const coexist = seg.coexistWith ?? [];
  if (coexist.length > 0)
    return `${seg.phase} + ${coexist.join(" + ")} (coexistence)`;
  return seg.phase;
}

export function PhaseStrip({
  segments,
  size = "md",
  emptyLabel = "No clear phase",
  orientation = "horizontal",
  className = "",
}: PhaseStripProps): JSX.Element {
  const vertical = orientation === "vertical";
  const indexed = segments
    .map((s) => s.phase)
    .filter((p): p is string => p !== null);
  const first = indexed.length > 0 ? indexed[0]! : null;
  const last = indexed.length > 0 ? indexed[indexed.length - 1]! : null;
  const distinct = new Set(indexed);

  return (
    <div className={cx(vertical && "h-full", className)} data-size={size} data-orientation={orientation}>
      <div className={cx("flex", vertical ? "h-full w-2 flex-col gap-[2px]" : sizeClass[size])}>
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
          const coexist = seg.coexistWith ?? [];
          return (
            <div
              key={i}
              data-testid="ps-seg"
              data-coexist-count={
                coexist.length > 0 ? String(coexist.length + 1) : undefined
              }
              aria-label={segLabel(seg)}
              title={segLabel(seg)}
              className="flex-1 rounded-[1.5px]"
              style={{ background: segBackground(seg) }}
            />
          );
        })}
      </div>
      {!vertical && (
      <div
        data-testid="ps-cap"
        className="mt-1.5 flex items-center gap-1.5 text-base text-ink-soft"
      >
        {first === null || last === null ? (
          <span className="font-semibold text-ink-soft">{emptyLabel}</span>
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
            <span className="text-ink-soft">throughout</span>
          </>
        )}
      </div>
      )}
    </div>
  );
}
