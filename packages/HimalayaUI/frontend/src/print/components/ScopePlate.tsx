import type { ReactNode } from "react";
import { Card, Kicker, PhaseStrip, Dot, Button, Field } from "../ui";
import type { PhaseSegment } from "../ui";
import { AutoGroup } from "./AutoGroup";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface ScopePlateProps {
  seriesName: string;
  /** Title is editable in-place; page-deferred (no input is built here). */
  onEditTitle?: () => void;
  /** AutoGroup body (with <strong> emphasis). */
  grouping: ReactNode;
  orderedBy: string;
  orderNote: ReactNode;
  onChangeOrder?: () => void;
  /** Ordering-variable options → makes the order-field a real dropdown. */
  orderOptions?: ReadonlyArray<{ value: string; label: ReactNode }>;
  onOrderSelect?: (value: string) => void;
  /** e.g. "6 samples · low to high". */
  count: string;
  onUndo?: () => void;
  undoLabel?: string;
  /** ScopeSampleRow×N (children-slotting). */
  rows: ReactNode;
  /** ScopeCandidateRow×N OR the empty node. */
  candidates: ReactNode;
  preview: PhaseSegment[];
  footState: { kind: "warn" | "ready"; text: string };
  footNote: ReactNode;
  buildDisabled?: boolean;
  /**
   * The scope→create chain is in flight: the foot button flips to the
   * progressive register ("Building…") with `aria-busy`, disabled (no
   * double-submit). The page derives this from the same stage/isPending
   * sources that gate the chain, so the label reverts on both terminal paths.
   */
  buildBusy?: boolean;
  onBuild?: () => void;
  /** Discard affordance, placed in the foot's action cluster immediately left
   *  of "Confirm & build" (SC-POLISH2). A ReactNode slot so the page owns the
   *  control AND its confirm-before-destroy state; ScopePlate only positions it.
   *  Omit it and no discard control renders. */
  discardSlot?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * ScopePlate — the whole series-scoping worksheet (mockup `.scope-plate`): a
 * centred, elevated review surface where the human confirms the machine's
 * grouping and ordering before building the series.
 *
 * Presentational: ScopePlate holds NO state. Every derivation — the preview
 * segments, the foot state line, the build gate, the count string, the undo
 * affordance — arrives as a prop computed by the consumer (the future Layer-4
 * page; simulated by ScopingAssembly). The row/candidate regions are ReactNode
 * slots so the page maps data → leaf rows. The rows wrapper strips the last
 * row's bottom hairline.
 */
export function ScopePlate({
  seriesName,
  onEditTitle,
  grouping,
  orderedBy,
  orderNote,
  onChangeOrder,
  orderOptions,
  onOrderSelect,
  count,
  onUndo,
  undoLabel,
  rows,
  candidates,
  preview,
  footState,
  footNote,
  buildDisabled,
  buildBusy,
  onBuild,
  discardSlot,
  className,
}: ScopePlateProps): JSX.Element {
  return (
    <Card elevated className={cx("w-full max-w-[760px] px-8 pt-7 pb-6", className)}>
      <Kicker tone="accent">New series</Kicker>

      {/* SC-POLISH2 (controls-don't-lie): the dotted-underline is the editable
          idiom. Render it ONLY when an edit handler is wired — otherwise the
          title is plain ink, so the affordance never promises an edit the page
          can't deliver. (The scoping worksheet's title is auto-derived and not
          yet editable; the builder edits its title through PlateTitle.) */}
      <h1 className="text-display text-ink mt-1.5">
        {onEditTitle ? (
          <span
            data-editable="true"
            className="border-b border-dotted border-hair-strong pb-px cursor-text"
            onClick={onEditTitle}
          >
            {seriesName}
          </span>
        ) : (
          <span>{seriesName}</span>
        )}
      </h1>

      <AutoGroup variant="summary" className="mt-4">
        {grouping}
      </AutoGroup>

      <Kicker as="h2" tone="soft" className="mt-5 mb-2">
        Ordered by
      </Kicker>
      <Field
        testId="order-field"
        srLabel="Ordered by"
        menuLabel="Ordered by"
        value={orderedBy}
        {...(orderOptions
          ? { options: orderOptions, ...(onOrderSelect ? { onSelect: onOrderSelect } : {}) }
          : onChangeOrder
            ? { onClick: onChangeOrder }
            : {})}
      />
      <div className="text-caption text-ink-soft mt-1.5">{orderNote}</div>

      <div className="flex items-baseline justify-between mt-6 mb-1">
        <Kicker as="h2" tone="soft">The series</Kicker>
        <div className="flex items-baseline gap-3.5">
          {onUndo ? (
            <button
              type="button"
              onClick={onUndo}
              className="text-caption font-semibold text-accent hover:underline"
              {...(undoLabel ? { title: undoLabel } : {})}
            >
              ↺ Undo last change
            </button>
          ) : null}
          <span className="text-caption font-mono text-ink-soft">{count}</span>
        </div>
      </div>

      <div className="[&>*:last-child]:border-b-0">{rows}</div>

      <div className="mt-5 pt-4 border-t border-hair">
        <Kicker as="h2" tone="soft" className="mb-2.5">
          Himalaya also found
        </Kicker>
        {candidates}
      </div>

      {/* Omitted entirely when nothing will be committed: a zero-segment
          PhaseStrip still paints an empty bar plus its "No clear phase"
          caption — a visible artifact previewing nothing. The foot's warn
          state already says "Keep at least one value to build". */}
      {preview.length > 0 ? (
        <div className="mt-5">
          <Kicker as="h2" tone="soft" className="mb-2">
            Preview · phase across the series
          </Kicker>
          <PhaseStrip segments={preview} size="sm" />
        </div>
      ) : null}

      <div className="mt-6 pt-4 border-t border-hair flex items-center justify-between gap-5">
        <div className="flex flex-col gap-1">
          <div
            className={cx(
              "flex items-center gap-2 text-meta font-semibold",
              footState.kind === "warn" ? "text-accent" : "text-ink",
            )}
          >
            <Dot tone={footState.kind === "warn" ? "accent" : "success"} aria-hidden />
            {footState.text}
          </div>
          <div className="text-caption text-ink-soft max-w-[42ch]">{footNote}</div>
        </div>
        <div className="flex items-center gap-2 flex-shrink-0">
          {discardSlot}
          <Button
            variant="solid"
            {...(buildDisabled || buildBusy ? { disabled: true } : {})}
            {...(buildBusy ? { "aria-busy": true } : {})}
            {...(onBuild && !buildBusy ? { onClick: onBuild } : {})}
          >
            {buildBusy ? "Building…" : "Confirm & build →"}
          </Button>
        </div>
      </div>
    </Card>
  );
}
