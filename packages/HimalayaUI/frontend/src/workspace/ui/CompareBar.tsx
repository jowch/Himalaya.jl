import { Kbd } from "./Kbd";
import type { CompareScale } from "../compare";
import { PLOT_VIEWS, type PlotView } from "../plotView";

export interface CompareChip {
  id: number;
  name: string;
  dash: string;
  scale: CompareScale | undefined;
}

interface Props {
  chips: readonly CompareChip[];
  subtractId: number | undefined;
  hidden: boolean;
  view: PlotView;
  onView: (view: PlotView) => void;
  onSubtract: (id: number) => void;
  onRemove: (id: number) => void;
  onAdd: () => void;
}

export function CompareBar({ chips, subtractId, hidden, view, onView, onSubtract, onRemove, onAdd }: Props): JSX.Element {
  return (
    <div className="flex flex-wrap items-center gap-2" data-testid="workspace-compare-bar">
      {chips.map((c) => {
        const sub = c.id === subtractId;
        return (
          <span
            key={c.id}
            data-compare-chip={c.id}
            className={
              "inline-flex items-center gap-1 rounded-sm border px-1.5 py-0.5 text-caption " +
              (sub ? "border-accent text-ink" : "border-hair text-ink-soft") +
              (hidden ? " opacity-50" : "")
            }
          >
            <svg width={18} height={6} aria-hidden="true">
              <line x1={0} x2={18} y1={3} y2={3} stroke="var(--color-ink-soft)" strokeWidth={1.4} strokeDasharray={c.dash} />
            </svg>
            <button
              type="button"
              aria-pressed={sub}
              title={sub ? "Stop subtracting" : "Subtract from the current trace"}
              onClick={() => onSubtract(c.id)}
              className="hover:text-ink"
            >
              {c.name}
            </button>
            {c.scale !== undefined && (
              <span
                className="text-ink-faint"
                data-scale-source={c.scale.source}
                title={c.scale.source === "monitor" ? "Scaled by the transmitted beam" : "No beam monitors, so matched at high q"}
              >
                ×{c.scale.scale.toPrecision(3)}
                {c.scale.source === "matched" && " matched"}
              </span>
            )}
            <button type="button" aria-label={`Remove ${c.name}`} onClick={() => onRemove(c.id)} className="text-ink-faint hover:text-ink">
              ×
            </button>
          </span>
        );
      })}
      <button type="button" onClick={onAdd} className="inline-flex items-center gap-1 text-caption text-ink-soft hover:text-ink">
        + Compare <Kbd>C</Kbd>
      </button>
      <span className="ml-auto inline-flex items-center gap-3 text-caption text-ink-faint">
        {chips.length > 0 && (
          <span className="inline-flex items-center gap-1">
            {hidden ? "Comparisons hidden" : "Comparisons shown"}
          </span>
        )}
        <span role="radiogroup" aria-label="Plot view" className="inline-flex items-center gap-2" data-testid="workspace-plot-view">
          {PLOT_VIEWS.map((o) => (
            <button
              key={o.view}
              type="button"
              role="radio"
              title={o.title}
              aria-checked={o.view === view}
              onClick={() => onView(o.view)}
              className={o.view === view ? "text-ink" : "hover:text-ink"}
            >
              {o.label}
            </button>
          ))}
        </span>
      </span>
    </div>
  );
}
