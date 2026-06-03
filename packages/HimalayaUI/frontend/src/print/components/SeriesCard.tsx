import type { ReactNode } from "react";
import { CardFigure } from "../waterfall/CardFigure";
import { PhaseStrip, NoticePill } from "../ui";
import type { PhaseSegment } from "../ui";
import type { WaterfallRow } from "../waterfall/waterfallModel";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface SeriesCardProps {
  /** Mini-waterfall rows, low→high; drives <CardFigure>. */
  rows: WaterfallRow[];
  /** Phase-strip segments, low→high; passed to <PhaseStrip>. */
  segments: PhaseSegment[];
  /** "Fig. 1" (saved) or "Recipe" (draft). */
  figLabel: string;
  title: string;
  /** Meta: "{n} samples". */
  sampleCount: number;
  /** Meta: "by {variable}". */
  variable: string;
  /** Footer-left: beamtime string or a cross-experiment node. */
  provenance: ReactNode;
  /** Footer-right: "2 days ago". */
  editedLabel: string;
  /** Footer-right author initials / name. */
  author: string;
  /** Optional kick-row pill. */
  notice?: { tone: "new"; count: number } | { tone: "draft" };
  /** Dashed border + dimmed figure (mockup .card.is-draft). */
  draft?: boolean;
  onClick?: () => void;
  /** PLACEMENT-ONLY, appended last. */
  className?: string;
}

/**
 * One card in the Series folio wall: a frozen mini-waterfall figure on top,
 * then a body with a kick-row (fig label + optional pill), a serif title, a
 * meta line, a phase strip, and a hairline footer (provenance + edit attribution).
 * A `draft` variant gets a dashed border and a dimmed figure.
 *
 * The closed-look/open-placement contract: appearance lives in the primitives
 * composed here; `className` is PLACEMENT-ONLY.
 */
export function SeriesCard({
  rows,
  segments,
  figLabel,
  title,
  sampleCount,
  variable,
  provenance,
  editedLabel,
  author,
  notice,
  draft = false,
  onClick,
  className,
}: SeriesCardProps): JSX.Element {
  return (
    <article
      data-testid="series-card"
      data-draft={draft ? "true" : "false"}
      onClick={onClick}
      className={cx(
        "bg-plate border border-hair rounded overflow-hidden",
        onClick && "cursor-pointer",
        draft && "border-dashed border-hair-strong",
        className,
      )}
    >
      {/* Figure region — the frozen mini-waterfall plate */}
      <div
        className={cx(
          "px-2 pt-2 border-b border-paper-sunk",
          draft && "opacity-60",
        )}
      >
        <CardFigure rows={rows} />
      </div>

      {/* Card body */}
      <div className="px-4 py-3">
        {/* Kick row: fig label + optional notice pill */}
        <div className="flex items-center justify-between mb-1">
          <span className="text-kicker text-kicker-accent">{figLabel}</span>
          {notice != null && (
            notice.tone === "new"
              ? <NoticePill tone="new">+{notice.count} new match</NoticePill>
              : <NoticePill tone="draft">Draft</NoticePill>
          )}
        </div>

        {/* Title — Newsreader serif, 19px via text-headline */}
        <h3 className="text-headline">{title}</h3>

        {/* Meta line */}
        <div className="flex items-center gap-1.5 mt-1 text-caption">
          <b className="text-ink-soft font-semibold">{sampleCount} samples</b>
          <span className="text-ink-faint" aria-hidden="true">·</span>
          <span className="text-ink-faint">by {variable}</span>
        </div>

        {/* Phase strip */}
        <PhaseStrip segments={segments} className="mt-3" />

        {/* Footer */}
        <div className="flex items-center justify-between mt-3 pt-2.5 border-t border-hair text-caption text-ink-faint">
          <span>{provenance}</span>
          <span>
            edited <b className="text-ink-soft font-semibold">{editedLabel}</b>
            {" · "}
            <b className="text-ink-soft font-semibold">{author}</b>
          </span>
        </div>
      </div>
    </article>
  );
}
