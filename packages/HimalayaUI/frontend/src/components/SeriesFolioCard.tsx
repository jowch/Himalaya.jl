import type { SeriesSummary } from "../api";
import { phaseColor } from "../phases";

interface SeriesFolioCardProps {
  series: SeriesSummary;
  /** Called with the series id when the card is activated. */
  onOpen: (id: number) => void;
}

const MAX_SWATCHES = 3;

/**
 * One masonry card in the series folio (#173 / I3.3). A lightweight "plate"
 * stand-in: a phase-swatch strip (proportional placeholder for the builder's
 * full mini-waterfall, deliberately deferred to #175), the title, author +
 * recency line, a member-count badge, draft-dashed styling, and a stale dot.
 *
 * Card height varies with member_count so the CSS-columns masonry reads as a
 * wall of distinct figures rather than a uniform grid (mockup intent).
 */
export function SeriesFolioCard({ series, onOpen }: SeriesFolioCardProps): JSX.Element {
  const isDraft = series.content_hash === "";
  const overflow = series.member_phase_count - series.member_phases.length;
  // Placeholder figure height scales with member_count (min 2 rows so a
  // zero-member draft still has a body). Tailwind can't take a dynamic
  // class; use an inline style for the computed height.
  const figRows = Math.max(series.member_count, 2);

  return (
    <button
      type="button"
      data-testid={`series-card-${series.id}`}
      data-member-count={series.member_count}
      data-draft={isDraft ? "true" : "false"}
      data-stale={series.has_stale_members ? "true" : "false"}
      onClick={() => onOpen(series.id)}
      className={[
        "mb-5 block w-full break-inside-avoid overflow-hidden rounded border text-left",
        "bg-plate transition-transform hover:-translate-y-0.5",
        isDraft ? "border-dashed border-hair-strong" : "border-hair",
      ].join(" ")}
    >
      {/* Placeholder "frozen plate" figure — proportional to member_count. */}
      <div
        data-testid={`series-card-${series.id}-fig`}
        className="flex flex-col justify-end gap-1 bg-paper-sunk p-3"
        style={{ minHeight: `${figRows * 18 + 24}px` }}
      >
        <div className="flex items-center gap-1" data-testid={`series-card-${series.id}-swatches`}>
          {series.member_phases.slice(0, MAX_SWATCHES).map((phase, i) => (
            <span
              key={`${phase}-${i}`}
              data-testid={`series-card-${series.id}-swatch-${i}`}
              title={phase}
              className="h-2.5 w-6 rounded-sm"
              style={{ background: phaseColor(phase) }}
            />
          ))}
          {overflow > 0 && (
            <span className="text-[10px] text-ink-faint">+{overflow} more</span>
          )}
        </div>
      </div>

      <div className="flex flex-col gap-1 p-3">
        <div className="flex items-start justify-between gap-2">
          <span className="text-sm font-semibold text-ink">
            {series.title.trim() === "" ? "Untitled series" : series.title}
          </span>
          {series.has_stale_members && (
            <span
              data-testid={`series-card-${series.id}-stale-dot`}
              title="Has stale members"
              className="mt-1 inline-block h-2 w-2 shrink-0 rounded-full bg-print-accent"
            />
          )}
        </div>
        <div className="flex items-center gap-2 text-xs text-ink-faint">
          <span>{series.author_username ?? "unknown"}</span>
          <span aria-hidden="true">·</span>
          <span data-testid={`series-card-${series.id}-count`}>
            {series.member_count} {series.member_count === 1 ? "member" : "members"}
          </span>
          {isDraft && (
            <>
              <span aria-hidden="true">·</span>
              <span className="uppercase tracking-wide">draft</span>
            </>
          )}
        </div>
      </div>
    </button>
  );
}
