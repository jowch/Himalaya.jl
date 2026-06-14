import type { ReactNode } from "react";
import { Button, IconButton, Kicker, Slider, Field, NoticePill } from "../ui";
import { AutoGroup } from "./AutoGroup";
import { RailSection } from "./RailSection";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface BuilderRailProps {
  grouping: ReactNode;
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
   * The "Edit" entry into draft state. Omit it (e.g. while a draft is already
   * live) and the affordance is NOT rendered — re-running the idempotent
   * ensureDraft is a no-op, and controls-don't-lie says we don't show it.
   */
  onAdjust?: () => void;
  /**
   * Discard the live draft (the "Cancel" verb). Present only in draft mode; it
   * pairs with "Save changes" as a proper primary/secondary BUTTON row, not a
   * link action buried under the recipe (BU-DRAFT-ACTIONS).
   */
  onCancel?: () => void;
  orderedBy: string;
  orderNote?: ReactNode;
  onChangeOrder?: () => void;
  /** Ordering-variable options → makes the field a real dropdown. */
  orderOptions?: ReadonlyArray<{ value: string; label: ReactNode }>;
  onOrderSelect?: (value: string) => void;
  offset: number;
  onOffsetChange: (v: number) => void;
  traces: ReactNode;
  /**
   * The traces slot is the editable recipe (draft mode), whose rows reorder via
   * drag + the ▲▼ buttons. Only then does the section label make the
   * "drag to reorder" promise — in read mode the slot is a static MemberList and
   * the instruction would lie (copy-doesn't-lie). Defaults to false (read mode).
   */
  reorderable?: boolean;
  onAddSample?: () => void;
  onCollapse?: () => void;
  /**
   * Foot caption. The default asserts WYSIWYG ("The plate above is the figure
   * as it will export…") — a read-state truth. A page in a state where that is
   * NOT true (e.g. mid-draft, when recipe edits haven't re-resolved the plate
   * yet) MUST override it with an honest variant (copy-doesn't-lie).
   */
  caption?: ReactNode;
  /** PLACEMENT-ONLY. Width is page-owned (the assembly sets it). */
  className?: string;
}

/**
 * BuilderRail — the series-builder "Compose" editing rail (mockup `.rail`).
 *
 * Presentational contract (Batch 7/9/10/11/12): holds NO local state. offset /
 * scale / the trace rows (the `traces` slot) / every handler are PROPS; the
 * SeriesBuilderAssembly page owns the real state and the figure↔rail link.
 *
 * LOAD-BEARING OMISSIONS ("controls don't lie"): the mockup's Representation
 * section (Waterfall/Heatmap) and the "Track reflections" toggle drive renderers
 * that are out-of-scope / deferred (HeatmapChart out-of-scope; TrackingLine
 * deferred — same call as Batch 9 SeriesPlate). Both are OMITTED; the Display
 * section keeps only the offset slider — the log/linear q-scale toggle lives on
 * the plate (a single contextual control, not a redundant rail+plate pair).
 * Figure export likewise lives on the PLATE head (`SeriesPlate`'s `actions`
 * slot, the ExportButton split button) — the rail-foot "Copy as PNG" was
 * removed for the same single-contextual-control reason.
 *
 * CONDITIONAL AFFORDANCES (same pattern as onAdjust): the "Collapse rail"
 * IconButton renders only when onCollapse is passed (the live builder page has
 * no collapse behavior yet), and the "+ Add sample" Button renders only when
 * onAddSample is passed (the page's real add path is the native select in the
 * traces slot); the foot row itself is conditional on onAddSample — when it is
 * absent, no foot row renders at all.
 */
export function BuilderRail({
  grouping,
  onConfirm,
  confirmBusy,
  confirmLabel,
  onAdjust,
  onCancel,
  orderedBy,
  orderNote,
  onChangeOrder,
  orderOptions,
  onOrderSelect,
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
        "flex flex-col gap-5 bg-paper-sunk border-l border-hair px-5 pt-4 pb-7 overflow-y-auto",
        className,
      )}
    >
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-2.5">
          <Kicker tone="soft">Compose</Kicker>
          {/* BU-MODESHIFT: draft (edit) mode is otherwise near-invisible against
              read mode (same bg/borders/geometry), so the editing state carries
              a standing "Unsaved draft" marker on the Compose header. `reorderable`
              is the draft signal (only a draft's traces reorder). */}
          {reorderable && <NoticePill tone="draft">Unsaved draft</NoticePill>}
        </div>
        {/* controls-don't-lie: no onCollapse means no collapse behavior exists,
            so the affordance is dropped; the Kicker sits flush-left on its own. */}
        {onCollapse && (
          <IconButton label="Collapse rail" tone="ghost" onClick={onCollapse}>
            &#8250;
          </IconButton>
        )}
      </div>

      {/* BU-AUTOGROUP-STALE: this series is user-owned and already saved, so the
          card carries NO "Auto-grouped" title. READ mode keeps the controls-don't-
          lie pair — a DISABLED "Save changes" (the affordance exists; nothing to
          save until you edit) + the muted "Edit" door. In DRAFT mode the card
          drops its link actions: the commit/discard decision renders as a proper
          primary/secondary BUTTON row below (BU-DRAFT-ACTIONS) instead of a link
          in a card plus a footnote-styled Cancel. */}
      <AutoGroup
        variant="compose"
        actions={
          onAdjust
            ? [
                // READ mode: Save is inert (no onClick) until a draft exists; Edit
                // is the live entry into draft state.
                { label: "Save changes" },
                { label: "Edit", muted: true, onClick: onAdjust },
              ]
            : []
        }
      >
        {grouping}
      </AutoGroup>

      {/* BU-DRAFT-ACTIONS: the live-draft commit/discard pair as real buttons.
          "Save changes" is the primary (accent); it states its busy reason while
          the Save→Commit chain runs and disables when the chain isn't ready. */}
      {reorderable && (
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
      )}

      <RailSection label="Ordering variable" {...(orderNote != null ? { note: orderNote } : {})}>
        <Field
          srLabel="Ordering variable"
          value={orderedBy}
          placeholder="Not set"
          {...(orderOptions
            ? { options: orderOptions, ...(onOrderSelect ? { onSelect: onOrderSelect } : {}) }
            : onChangeOrder
              ? { onClick: onChangeOrder }
              : {})}
        />
      </RailSection>

      {/* DISPLAY holds the Trace-offset slider only. The log/linear-q scale
          toggle lives on the PLATE (contextual to the figure that exports) —
          a single q-scale control, not a redundant pair. */}
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

      {/* controls-don't-lie: the live page adds samples through the native
          select in the traces slot, so the foot row only renders when a rail
          add path (onAddSample) actually exists. */}
      {onAddSample && (
        <div className="flex gap-2">
          <Button variant="outline" className="flex-1" onClick={onAddSample}>
            + Add sample
          </Button>
        </div>
      )}

      <div className="text-caption text-ink-soft leading-relaxed">
        {caption ??
          "The plate above is the figure as it will export. What you compose is what you publish."}
      </div>
    </aside>
  );
}
