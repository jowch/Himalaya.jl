import type { ReactNode } from "react";
import { RepresentationToggle, type Representation } from "./RepresentationToggle";

interface SeriesBuilderRailProps {
  collapsed: boolean;
  onToggleCollapsed: () => void;
  representation: Representation;
  onRepresentationChange: (next: Representation) => void;
  orderingVariable: string | null;
  /** Parent injects <FigureExportControls/> (it owns the export spec thunk). */
  exportControls: ReactNode;
}

/**
 * SeriesBuilderRail — the quiet editing margin of the series builder (#175).
 * Read-only meta + the representation toggle + injected export controls. The
 * rail collapses to full-bleed (the mockup's only "third representation" is
 * this rail folding away, not a distinct plot mode). Recipe-edit controls
 * (add sample, reorder commit, ordering-variable picker) are I3.5b — this
 * surface only displays.
 */
export function SeriesBuilderRail({
  collapsed, onToggleCollapsed, representation, onRepresentationChange,
  orderingVariable, exportControls,
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
      className="flex w-[336px] shrink-0 flex-col gap-5 overflow-y-auto border-l border-hair p-4"
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

      <section className="flex flex-col gap-1">
        <div className="text-xs font-semibold text-ink-faint">Ordering variable</div>
        <div className="text-sm text-ink">{orderingVariable ?? "—"}</div>
      </section>

      <section className="flex flex-col gap-1.5">
        <div className="text-xs font-semibold text-ink-faint">Representation</div>
        <RepresentationToggle value={representation} onChange={onRepresentationChange} />
      </section>

      <section className="flex flex-col gap-1.5" data-testid="rail-export">
        {exportControls}
      </section>
    </aside>
  );
}
