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
  /** Autogroup card (B-I): how many samples Himalaya read as this series. */
  sampleCount: number;
  onConfirmSeries: () => void;
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
  sampleCount, onConfirmSeries, onAdjustSeries, exportControls, editControls,
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
        <span className="text-xs font-semibold uppercase tracking-wide text-ink">Compose</span>
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
        <section className="flex flex-col gap-1.5 [&_input]:bg-plate" data-testid="rail-edit">
          {editControls}
        </section>
      ) : (
        <>
          <AutogroupCard
            sampleCount={sampleCount}
            orderingVariable={orderingVariable}
            onConfirm={onConfirmSeries}
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
      </section>

      <section className="flex flex-col gap-1.5" data-testid="rail-export">
        {exportControls}
      </section>
    </aside>
  );
}
