import type { ReactNode } from "react";
import { Button, IconButton, Kicker, Slider, NoticePill } from "../ui";
import { RailSection } from "./RailSection";
import { cx } from "../../lib/cx";


export interface BuilderRailProps {
  onConfirm?: () => void;
  /**
   * The Save→Commit chain is in flight: the "Save changes" action flips to the
   * progressive register ("Saving…") with `aria-busy`, still disabled — the
   * control states WHY it is inert instead of silently sitting dead. The page
   * derives this from the same stage/isPending sources that gate the chain, so
   * the label reverts on both terminal paths.
   */
  confirmBusy?: boolean;
  /**
   * BU-PROGRESS: the step-named busy label while the Save→Commit chain runs
   * ("Saving order…" / "Publishing the figure…" / "Confirming…"). When busy and
   * present it replaces the generic "Saving…" so the action shows visible
   * progress through the multi-step chain. Falls back to "Saving…".
   */
  confirmLabel?: string;
  /**
   * The "Edit" entry into edit mode. Read mode shows ONLY this; while a draft is
   * live it is withheld and the Save changes / Cancel pair shows instead.
   */
  onAdjust?: () => void;
  /**
   * Discard the live draft (the "Cancel" verb). Present only in edit mode; it
   * pairs with "Save changes" in the footer's primary/secondary button row.
   */
  onCancel?: () => void;
  offset: number;
  onOffsetChange: (v: number) => void;
  traces: ReactNode;
  /**
   * Edit-mode flag (a live draft). It drives the "Unsaved draft" marker, the
   * "drag to reorder" Traces label (the read-mode MemberList cannot be dragged,
   * so the promise would lie), and which footer action pair shows.
   */
  reorderable?: boolean;
  onAddSample?: () => void;
  onCollapse?: () => void;
  /**
   * Footer caption. The default asserts WYSIWYG ("The plate above is the figure
   * as it will export…") — a read-state truth. Mid-draft (when recipe edits have
   * not re-resolved the plate yet) the page overrides it with an honest variant.
   */
  caption?: ReactNode;
  /** PLACEMENT-ONLY. Width is page-owned (the assembly sets it). */
  className?: string;
}

/**
 * BuilderRail — the series-builder "Compose" panel.
 *
 * Two-zone layout: a scrolling content column (Display + Traces) over a pinned
 * footer that holds the WYSIWYG caption and the mode-switch actions. Read mode
 * is genuinely read-only — the only affordance is the footer "Edit" door;
 * editing the title, reordering, and adding samples all live behind it (a real
 * mode toggle, not a lazy auto-draft). A live draft swaps the footer to the
 * Save changes / Cancel pair and unlocks the editable recipe in the traces slot.
 *
 * Presentational contract: holds NO local state. offset / the trace rows (the
 * `traces` slot) / every handler are PROPS; SeriesBuilderPage owns the state and
 * the figure↔rail link.
 *
 * LOAD-BEARING OMISSIONS ("controls don't lie"): the mockup's Representation
 * section (Waterfall/Heatmap) and the "Track reflections" toggle drive renderers
 * that are out-of-scope / deferred, so they are omitted. The Display section
 * keeps only the offset slider — the log/linear q-scale toggle and figure export
 * live on the PLATE (single contextual controls, not redundant rail+plate pairs).
 * The Ordering-variable readout was dropped: it is set during scoping and is not
 * editable here, so a static field only added chrome.
 */
export function BuilderRail({
  onConfirm,
  confirmBusy,
  confirmLabel,
  onAdjust,
  onCancel,
  offset,
  onOffsetChange,
  traces,
  reorderable = false,
  onAddSample,
  onCollapse,
  caption,
  className,
}: BuilderRailProps): JSX.Element {
  return (
    <aside
      data-testid="builder-rail"
      className={cx(
        "flex flex-col bg-paper-sunk border-l border-hair",
        className,
      )}
    >
      {/* Scrolling content: header + the view/edit sections. */}
      <div className="flex flex-col gap-5 px-5 pt-4 pb-6 overflow-y-auto grow">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2.5">
            <Kicker tone="soft">Compose</Kicker>
            {/* BU-MODESHIFT: edit mode is otherwise near-invisible against read
                mode (same surface/geometry), so it carries a standing "Unsaved
                draft" marker. `reorderable` is the edit-mode signal. */}
            {reorderable && <NoticePill tone="draft">Unsaved draft</NoticePill>}
          </div>
          {/* controls-don't-lie: no onCollapse means no collapse behavior. */}
          {onCollapse && (
            <IconButton label="Collapse rail" tone="ghost" onClick={onCollapse}>
              &#8250;
            </IconButton>
          )}
        </div>

        {/* DISPLAY holds the Trace-offset slider only — a view control, live in
            both modes (it never edits the series). */}
        <RailSection label="Display">
          <Slider
            label="Trace offset"
            valueDisplay={`${offset.toFixed(2)}×`}
            value={offset}
            min={0.4}
            max={1.4}
            step={0.05}
            onChange={onOffsetChange}
          />
        </RailSection>

        {/* copy-doesn't-lie: only the draft (reorderable) traces slot can be
            dragged, so only then does the label make the reorder promise. */}
        <RailSection label={reorderable ? "Traces · drag to reorder" : "Traces"}>
          <div className="flex flex-col gap-0.5">{traces}</div>
        </RailSection>

        {/* controls-don't-lie: the live page adds samples through the picker in
            the traces slot, so this rail-foot add path only renders when wired. */}
        {onAddSample && (
          <Button variant="outline" className="w-full" onClick={onAddSample}>
            + Add sample
          </Button>
        )}
      </div>

      {/* Pinned footer: the WYSIWYG caption and the mode-switch actions. A live
          draft (reorderable) shows Save changes (accent) + Cancel; read mode shows
          the single Edit door — the only editing entry point. */}
      <div className="border-t border-hair px-5 py-4 flex flex-col gap-3">
        <p className="text-caption text-ink-soft leading-relaxed">
          {caption ??
            "The plate above is the figure as it will export. What you compose is what you publish."}
        </p>
        {reorderable ? (
          <div className="flex gap-2">
            <Button
              variant="accent"
              className="flex-1"
              data-testid="builder-save"
              disabled={onConfirm === undefined}
              {...(onConfirm ? { onClick: onConfirm } : {})}
              {...(confirmBusy ? { "aria-busy": true } : {})}
            >
              {confirmBusy ? (confirmLabel ?? "Saving…") : "Save changes"}
            </Button>
            {onCancel && (
              <Button variant="outline" data-testid="builder-cancel" onClick={onCancel}>
                Cancel
              </Button>
            )}
          </div>
        ) : onAdjust ? (
          <Button
            variant="outline"
            className="w-full"
            data-testid="builder-edit"
            onClick={onAdjust}
          >
            Edit
          </Button>
        ) : null}
      </div>
    </aside>
  );
}
