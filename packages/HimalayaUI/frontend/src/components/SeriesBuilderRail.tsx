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
  /**
   * I3.5b — recipe-edit controls (add/remove/reorder sample, ordering-var,
   * order-rule, Save/Commit/Cancel). Injected by the page only in edit mode;
   * when absent the rail is the read-only display from I3.5a (the static
   * ordering-variable line below renders instead).
   */
  editControls?: ReactNode;
}

/**
 * SeriesBuilderRail — the editing margin of the series builder. In read mode
 * (#175 / I3.5a) it shows static meta + the representation toggle + injected
 * export controls. In edit mode (I3.5b) the page injects `editControls` (the
 * recipe list + add-sample picker + ordering-var/order-rule + Save/Commit).
 * The rail collapses to full-bleed (the mockup's only "third representation").
 */
export function SeriesBuilderRail({
  collapsed, onToggleCollapsed, representation, onRepresentationChange,
  orderingVariable, exportControls, editControls,
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

      {editControls !== undefined ? (
        <section className="flex flex-col gap-1.5" data-testid="rail-edit">
          {editControls}
        </section>
      ) : (
        <section className="flex flex-col gap-1">
          <div className="text-xs font-semibold text-ink-faint">Ordering variable</div>
          <div className="text-sm text-ink">{orderingVariable ?? "—"}</div>
        </section>
      )}

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
