import type { SeriesSummary } from "../api";
import { phaseColor } from "../phases";
import { useSeries } from "../queries";
import { SeriesMiniWaterfall } from "./SeriesMiniWaterfall";
import { PhaseStrip } from "./ui/PhaseStrip";
import { Kicker } from "./ui/Kicker";
import { buildPhaseStrip } from "../lib/series/folioFigure";

interface SeriesFolioCardProps {
  series: SeriesSummary;
  /** Called with the series id when the card is activated. */
  onOpen: (id: number) => void;
  /** 1-based figure ordinal among committed cards (shown in the kicker). */
  figNumber?: number;
  /**
   * Count of corpus samples that match this series' recipe but aren't members
   * yet (the "+N new match" pill, finding F-C). Not on the listing yet — the
   * pill is wired and defaults to 0; the caller supplies a value once a
   * recipe-match source exists.
   */
  newMatches?: number;
  /** Injectable clock for deterministic edited-timestamp tests. */
  now?: Date;
}

const MAX_SWATCHES = 3;

/** Relative-time formatting matching the mockup's `edited` ladder. */
function formatEdited(iso: string | null, now: Date): string {
  if (iso === null) return "recently";
  // SQLite `datetime()` strings are space-separated UTC; normalise to ISO.
  const then = new Date(iso.replace(" ", "T") + "Z");
  if (Number.isNaN(then.getTime())) return "recently";
  const days = Math.floor((now.getTime() - then.getTime()) / 86_400_000);
  if (days <= 0) return "just now";
  if (days === 1) return "yesterday";
  if (days < 7) return `${days} days ago`;
  if (days < 14) return "1 week ago";
  if (days < 21) return "2 weeks ago";
  const weeks = Math.floor(days / 7);
  return `${weeks} weeks ago`;
}

/**
 * One masonry card in the series folio (#173 / I3.3; R6 #229). A real plate:
 * a LIVE per-series mini-waterfall + a per-sample phase strip (read from the
 * series detail's ordered member snapshots), a serif title, a recipe/draft
 * pill, an ordering-variable meta line, and a provenance footer.
 *
 * The card pulls its full detail (`useSeries`) for the ordered per-sample data
 * the corpus listing doesn't carry. While that's loading — or for a draft with
 * no members yet — it falls back to the lightweight phase-swatch strip derived
 * from the listing's top-3 phases, so the card is never blank.
 *
 * Card height varies with member_count (the figure is taller for longer
 * series), so the CSS-columns masonry reads as a wall of distinct figures.
 */
export function SeriesFolioCard({
  series,
  onOpen,
  figNumber,
  newMatches = 0,
  now = new Date(),
}: SeriesFolioCardProps): JSX.Element {
  const isDraft = series.content_hash === "";
  const detail = useSeries(series.id);
  const members = detail.data?.members ?? [];
  const hasMiniature = members.length > 0;

  const overflow = series.member_phase_count - series.member_phases.length;
  const figRows = Math.max(series.member_count, 2);

  const orderingVariable = detail.data?.ordering_variable ?? null;

  const pill = isDraft
    ? { kind: "draft", text: "Draft" }
    : newMatches > 0
      ? { kind: "new", text: `+${newMatches} new match` }
      : null;

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
      {/* The frozen plate — a live miniature once detail loads, else a swatch
          strip stand-in (cold/loading or a zero-member draft). */}
      <div
        data-testid={`series-card-${series.id}-fig`}
        className={[
          "bg-paper-sunk p-2",
          isDraft ? "opacity-[0.62]" : "",
        ].join(" ")}
        style={hasMiniature ? undefined : { minHeight: `${figRows * 18 + 24}px` }}
      >
        {hasMiniature ? (
          <SeriesMiniWaterfall members={members} />
        ) : (
          <div
            className="flex h-full items-end gap-1"
            data-testid={`series-card-${series.id}-swatches`}
          >
            {series.member_phases.slice(0, MAX_SWATCHES).map((phase, i) => (
              <span
                key={`${phase}-${i}`}
                data-testid={`series-card-${series.id}-swatch-${i}`}
                title={phase}
                className="h-2.5 w-6 rounded-sm"
                style={{ background: phaseColor(phase) }}
              />
            ))}
            {overflow > 0 && <span className="text-[10px] text-ink-faint">+{overflow} more</span>}
          </div>
        )}
      </div>

      <div className="flex flex-col gap-1 p-3">
        {/* kicker + pill row */}
        <div className="flex items-center justify-between gap-2">
          <Kicker as="span" tone="accent" data-testid={`series-card-${series.id}-kicker`}>
            {isDraft ? "Recipe" : `Fig. ${figNumber ?? series.id}`}
          </Kicker>
          {pill && (
            <span
              data-testid={`series-card-${series.id}-pill`}
              className={[
                "rounded-full px-2 py-0.5 text-[10px] font-bold tracking-[0.03em]",
                pill.kind === "draft"
                  ? "border border-dashed border-hair-strong bg-paper-sunk text-ink-faint"
                  : "text-print-accent",
              ].join(" ")}
              style={
                pill.kind === "new"
                  ? { background: "color-mix(in oklab, var(--color-print-accent) 14%, transparent)" }
                  : undefined
              }
            >
              {pill.text}
            </span>
          )}
        </div>

        {/* serif title (X-1) */}
        <div className="flex items-start justify-between gap-2">
          <span className="text-headline leading-tight text-ink">
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

        {/* meta line */}
        <div
          data-testid={`series-card-${series.id}-meta`}
          className="flex items-center gap-1.5 text-[11.5px] text-ink-faint"
        >
          <b className="font-semibold text-ink-soft">
            {series.member_count} {series.member_count === 1 ? "sample" : "samples"}
          </b>
          {orderingVariable && (
            <>
              <span aria-hidden="true">·</span>
              <span>by {orderingVariable}</span>
            </>
          )}
        </div>

        {/* live per-sample phase strip + caption (F-B); only with detail */}
        {hasMiniature && <PhaseStrip segments={buildPhaseStrip(members).segments} className="mt-3" />}

        {/* footer rule + provenance + edited timestamp (F-H) */}
        <div
          data-testid={`series-card-${series.id}-foot`}
          className="mt-3 flex items-center justify-between border-t border-hair pt-2.5 text-[10.5px] text-ink-faint"
        >
          <span>
            {series.member_count} {series.member_count === 1 ? "member" : "members"}
            {orderingVariable ? ` · ${orderingVariable}` : ""}
          </span>
          <span>
            edited <b className="font-semibold text-ink-soft">{formatEdited(series.updated_at ?? series.last_event_at, now)}</b>
            {series.author_username ? ` · ${series.author_username}` : ""}
          </span>
        </div>
      </div>
    </button>
  );
}
