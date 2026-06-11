import { Fragment } from "react";
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
 * Caption truthfulness (SC-PREVCLAIM): "throughout" may only claim what every
 * cell shows. A single-phase strip with unindexed gaps (null / form_factor /
 * null-state cells) carries the coverage fraction instead ("Pn3m in 1 of 3"),
 * and single-phase coexistence is named, not erased ("Im3m + Lamellar in
 * 1 of 3"; transition captions keep their bare "A → B" shape — partners
 * there live in the segment aria/titles).
 * "Throughout" survives a coexistence strip only when every cell shows the
 * identical story (same phase, same partners) — otherwise the fraction form is
 * used even at full coverage ("in 3 of 3"), because the partner is present
 * somewhere but not everywhere and either bare claim would overstate.
 * Transition captions keep their "A → B" shape and append the same fraction
 * when coverage is partial.
 *
 * Per-segment `role="img"` + `aria-label`/`title` carry the phase name as the
 * accessible second channel (colour is never the sole signal — see phases.ts);
 * the role follows the house Semantic-Dot pattern (Dot.tsx) so the label is
 * reliably exposed. The caption's visual transition arrow is decorative
 * (aria-hidden) with an sr-only "to" carrying the relation for AT. The visual
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
  // Coverage honesty (SC-PREVCLAIM): the caption may not say more than the
  // cells show. `partial` → the fraction form replaces / accompanies the claim.
  const total = segments.length;
  const partial = indexed.length < total;
  // Distinct coexisting partners beyond each cell's dominant phase, in
  // appearance order — a coexistence strip never captions the dominant alone.
  const partners: string[] = [];
  for (const seg of segments) {
    if (seg.phase === null) continue;
    for (const p of seg.coexistWith ?? []) {
      if (p !== seg.phase && !partners.includes(p)) partners.push(p);
    }
  }
  // "Throughout" with partners requires every cell to tell the identical
  // story (same phase, same partner SET — sorted, so partner order cannot
  // fake non-uniformity) — see the header note.
  const segStory = (s: PhaseSegment): string =>
    `${s.phase}|${[...(s.coexistWith ?? [])].sort().join("+")}`;
  const uniform =
    total > 0 && segments.every((s) => segStory(s) === segStory(segments[0]!));
  const fraction = (
    <span className="text-ink-soft">
      in {indexed.length} of {total}
    </span>
  );

  return (
    <div className={cx(vertical && "h-full", className)} data-size={size} data-orientation={orientation}>
      <div className={cx("flex", vertical ? "h-full w-2 flex-col gap-[2px]" : sizeClass[size])}>
        {segments.map((seg, i) => {
          // Form-factor → hollow dashed cell; null → a faint distinct cell.
          if (seg.state === "form_factor") {
            return (
              <div
                key={i}
                role="img"
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
                role="img"
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
              role="img"
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
            <span className="sr-only">to</span>
            <span className="font-semibold" style={{ color: phaseColor(last) }}>
              {last}
            </span>
            {partial && fraction}
          </>
        ) : (
          <>
            <span className="font-semibold" style={{ color: phaseColor(first) }}>
              {first}
            </span>
            {/* Coexistence marker: "+" is decorative, sr says "with" (the
                same pattern as the transition arrow's sr-only "to"). */}
            {partners.map((p) => (
              <Fragment key={p}>
                <span className="text-ink-faint" aria-hidden="true">
                  +
                </span>
                <span className="sr-only">with</span>
                <span className="font-semibold" style={{ color: phaseColor(p) }}>
                  {p}
                </span>
              </Fragment>
            ))}
            {!partial && (partners.length === 0 || uniform) ? (
              <span className="text-ink-soft">throughout</span>
            ) : (
              fraction
            )}
          </>
        )}
      </div>
      )}
    </div>
  );
}
