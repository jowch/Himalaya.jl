import type { ReactNode } from "react";
import { RepresentationToggle, type Representation } from "./RepresentationToggle";
import { OffsetSlider } from "./OffsetSlider";
import { ScaleToggle, type ScaleMode } from "./ScaleToggle";
import { AutogroupCard } from "./AutogroupCard";

interface SeriesBuilderRailProps {
  collapsed: boolean;
  onToggleCollapsed: () => void;
  representation: Representation;
  onRepresentationChange: (next: Representation) => void;
  orderingVariable: string | null;
  /** Trace-offset slider (B-F). 0.4..1.4. */
  offset: number;
  onOffsetChange: (next: number) => void;
  /** log/linear q-axis scale (B-F). */
  scaleMode: ScaleMode;
  onScaleModeChange: (next: ScaleMode) => void;
  /** Cross-trace peak-tracking annotation (#208). Mockup `series-builder.html:469-476`. */
  trackOn: boolean;
  onTrackOnChange: (next: boolean) => void;
  /** Autogroup card (B-I): how many samples Himalaya read as this series. */
  sampleCount: number;
  onAdjustSeries: () => void;
  /** Parent injects <FigureExportControls/> (it owns the export spec thunk). */
  exportControls: ReactNode;
  /**
   * I3.5b — recipe-edit controls (add/remove/reorder sample, ordering-var,
   * order-rule, Save/Commit/Cancel). Injected only in edit mode; when present
   * the autogroup card + the static ordering-variable line are replaced by the
   * recipe editor. The compose controls (offset slider + scale toggle) stay in
   * both modes — they shape the figure regardless of whether the recipe is
   * being edited.
   */
  editControls?: ReactNode;
}

/**
 * SeriesBuilderRail — the editing margin of the series builder. Read mode
 * (#175 / I3.5a) shows the autogroup card (B-I), the static ordering-variable
 * line, the representation toggle, the compose controls (offset slider +
 * scale toggle, B-F), and the injected export controls. Edit mode (I3.5b)
 * swaps the autogroup + ordering line for the recipe editor.
 *
 * The rail collapses to full-bleed (the mockup's only "third representation");
 * the floating OffsetDock (rendered by the page) keeps offset reachable while
 * collapsed (B-G).
 *
 * B-K/B-L: edit inputs use `bg-plate`; the rail itself is the recessed
 * `bg-paper-sunk` margin against the bright plate (mockup
 * `.rail{background:var(--paper-sunk)}`).
 */
export function SeriesBuilderRail({
  collapsed, onToggleCollapsed, representation, onRepresentationChange,
  orderingVariable, offset, onOffsetChange, scaleMode, onScaleModeChange,
  trackOn, onTrackOnChange,
  sampleCount, onAdjustSeries, exportControls, editControls,
}: SeriesBuilderRailProps): JSX.Element {
  if (collapsed) {
    return (
      <button
        type="button"
        data-testid="rail-restore"
        onClick={onToggleCollapsed}
        title="Show the editing rail"
        className="flex items-center gap-1.5 border-l border-hair px-2 text-xs text-ink-faint"
      >
        <span aria-hidden="true">‹</span> Compose
      </button>
    );
  }
  return (
    <aside
      data-testid="series-builder-rail"
      className="flex w-[336px] shrink-0 flex-col gap-5 overflow-y-auto border-l border-hair bg-paper-sunk p-4"
    >
      <div className="flex items-center justify-between">
        {/*
          R3-Y03: the rail recedes behind the figure plate (DESIGN.md §4
          Flat-Except-the-Plate Rule). Ink-faint + the kicker letter-spacing,
          no semibold — a quiet section label, not a title competing with the
          plate's terracotta kicker.
        */}
        <span
          data-testid="rail-compose-header"
          data-recede="true"
          className="text-xs uppercase tracking-[0.14em] text-ink-faint"
        >Compose</span>
        <button
          type="button"
          data-testid="rail-collapse-toggle"
          onClick={onToggleCollapsed}
          title="Collapse the rail — full-bleed"
          className="rounded px-1.5 text-ink-faint hover:text-ink"
        >
          ›
        </button>
      </div>

      {editControls !== undefined ? (
        <section
          className="flex flex-col gap-1.5 rail-edit-inputs"
          data-testid="rail-edit"
          data-rail-edit-inputs=""
        >
          {/*
            R3-Y09: the plate-fill applies to text inputs only via the scoped
            `.rail-edit-inputs` rule in styles.css — the old `[&_input]`
            wildcard also caught injected `<input type="range">` slider thumbs.
          */}
          {editControls}
        </section>
      ) : (
        <>
          <AutogroupCard
            sampleCount={sampleCount}
            orderingVariable={orderingVariable}
            onAdjust={onAdjustSeries}
          />
          <section className="flex flex-col gap-1">
            <div className="text-xs font-semibold text-ink-faint">Ordering variable</div>
            <div className="text-sm text-ink">{orderingVariable ?? "—"}</div>
          </section>
        </>
      )}

      <section className="flex flex-col gap-1.5">
        <div className="text-xs font-semibold text-ink-faint">Representation</div>
        <RepresentationToggle value={representation} onChange={onRepresentationChange} />
      </section>

      <section className="flex flex-col gap-2.5" data-testid="rail-display">
        <div className="text-xs font-semibold text-ink-faint">Display</div>
        <OffsetSlider value={offset} onChange={onOffsetChange} />
        <ScaleToggle value={scaleMode} onChange={onScaleModeChange} />
        {/*
          Cross-trace peak-tracking (#208). The mockup gates this behind a
          checkbox: it reads best when the user already has a phase call and
          wants to see migration. Default off so it never crowds an unindexed
          stack. Waterfall-only in the mockup; we let it render in heatmap
          too because the y-band envelope is the same, but hide the row when
          there's nothing to connect (the layer self-empties).
        */}
        <label
          data-testid="track-toggle"
          className="flex items-center gap-2 text-xs text-ink-soft"
        >
          <input
            type="checkbox"
            data-testid="track-toggle-input"
            data-accent="print"
            checked={trackOn}
            onChange={(e) => onTrackOnChange(e.target.checked)}
            className="rounded border-hair accent-print-accent"
          />
          Track reflections
        </label>
      </section>

      <section className="flex flex-col gap-1.5" data-testid="rail-export">
        {exportControls}
      </section>
    </aside>
  );
}
