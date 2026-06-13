import type { ReactNode } from "react";
import { CardFigure } from "../waterfall/CardFigure";
import { Card, PhaseStrip, NoticePill, Button } from "../ui";
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
  /** Controls what the figure region renders.
   *  - `"ready"` (default): CardFigure + PhaseStrip (normal path).
   *  - `"pending"`: skeleton bone; strip suppressed (no false "No clear phase").
   *  - `"error"`: honest note; strip suppressed. */
  figureState?: "ready" | "pending" | "error";
  /** Optional retry handler for the `"error"` figure tile (F4). When provided,
   *  the error tile shows a "Try again" control that re-triggers this card's
   *  figure fetch — no full-page reload. Omitted ⇒ no control (the tile stays
   *  the honest note only; existing consumers stay byte-identical). */
  onRetryFigure?: () => void;
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
  figureState = "ready",
  onRetryFigure,
  onClick,
  className,
}: SeriesCardProps): JSX.Element {
  // The plate look comes from the Card primitive (the single source of the
  // Plate Lift — DESIGN.md §4: folio cards ARE lifted plates). Saved cards are
  // `elevated`; drafts take Card's dashed flat variant (no shadow — a recipe,
  // not yet a print). `interactive` adds the quiet house hover affordance
  // (hairline firming, no motion) when the card is a door.
  return (
    <Card
      as="article"
      elevated={!draft}
      draft={draft}
      interactive={onClick !== undefined}
      data-testid="series-card"
      data-draft={draft ? "true" : "false"}
      data-figure-state={figureState}
      {...(onClick ? { onClick } : {})}
      className={cx("overflow-hidden", className)}
    >
      {/* Figure region — the frozen mini-waterfall plate */}
      <div
        className={cx(
          "px-2 pt-2 border-b border-paper-sunk",
          draft && "opacity-60",
        )}
      >
        {figureState === "ready" && <CardFigure rows={rows} />}
        {figureState === "pending" && (
          <div
            data-testid="card-figure-pending"
            className="h-28 rounded-md bg-paper-sunk"
            aria-hidden="true"
          />
        )}
        {figureState === "error" && (
          <div
            data-testid="card-figure-error"
            className="h-28 flex flex-col items-center justify-center gap-1.5 text-caption text-ink-soft"
          >
            <span>Couldn&apos;t load this figure</span>
            {onRetryFigure !== undefined && (
              // stopPropagation: the retry must not also fire the card's
              // whole-card navigation (the article onClick door).
              <Button
                variant="outline"
                onClick={(e) => {
                  e.stopPropagation();
                  onRetryFigure();
                }}
              >
                Try again
              </Button>
            )}
          </div>
        )}
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

        {/* Title — Newsreader serif, 19px via text-headline.
            When the card is interactive (onClick present), the title hosts a
            real <button> so keyboard/SR users have a focusable primary action
            (Tab → Enter/Space) — the whole-card <article onClick> stays a
            mouse-only enhancement (FOL-KBD, WCAG 2.1.1 / 2.4.7). stopPropagation
            prevents the title click also bubbling to the article (double-nav).
            With no onClick, the title is plain text — no dead control. */}
        <h2 className="text-headline">
          {onClick ? (
            <button
              type="button"
              onClick={(e) => {
                e.stopPropagation();
                onClick();
              }}
              className="text-left focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent"
            >
              {title}
            </button>
          ) : (
            title
          )}
        </h2>

        {/* Meta line — the "· by {variable}" ordering clause is dropped when the
            series has no ordering variable (e.g. a form-factor series), so the
            line never reads "… · by" with a dangling preposition. */}
        <div className="flex items-center gap-1.5 mt-1 text-caption">
          <b className="text-ink-soft font-semibold">{sampleCount} samples</b>
          {variable.trim() !== "" && (
            <>
              <span className="text-ink-faint" aria-hidden="true">·</span>
              <span className="text-ink-soft">by {variable}</span>
            </>
          )}
        </div>

        {/* Phase strip — suppressed when figure is not ready to avoid the false
            "No clear phase" affirmative claim (FOL-HONEST-DERIVED). */}
        {figureState === "ready" && <PhaseStrip segments={segments} className="mt-3" />}

        {/* Footer */}
        <div className="flex items-center justify-between mt-3 pt-2.5 border-t border-hair text-caption text-ink-soft">
          <span>{provenance}</span>
          <span>
            edited <b className="text-ink-soft font-semibold">{editedLabel}</b>
            {/* FOL-RESCORE3: drop the "· author" clause when there is no author,
                so the line never trails a dangling separator. */}
            {author ? (
              <>
                {" · "}
                <b className="text-ink-soft font-semibold">{author}</b>
              </>
            ) : null}
          </span>
        </div>
      </div>
    </Card>
  );
}
