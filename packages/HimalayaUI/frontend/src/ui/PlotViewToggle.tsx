import { Tooltip } from "./Tooltip";
import type { PlotView, PlotViewOption } from "../workspace/plotView";

interface Props {
  options: readonly PlotViewOption[];
  view: PlotView;
  onChange: (view: PlotView) => void;
}

export function PlotViewToggle({ options, view, onChange }: Props): JSX.Element {
  return (
    <div
      role="group"
      aria-label="Plot view"
      data-testid="workspace-plot-view"
      className="flex h-5 shrink-0 items-center rounded-sm border border-hair bg-paper/90 text-caption leading-none"
    >
      {options.map((o) => (
        <Tooltip key={o.view} content={o.short === o.name ? "Plot view" : `Plot view: ${o.name}`} shortcut="V">
          <button
            type="button"
            aria-pressed={o.view === view}
            {...(o.short === o.name ? {} : { "aria-label": o.name })}
            onClick={() => onChange(o.view)}
            className={"h-full whitespace-nowrap rounded-sm px-1.5 " + (o.view === view ? "bg-plate text-ink shadow-sm" : "text-ink-soft hover:text-ink")}
          >
            {o.short}
          </button>
        </Tooltip>
      ))}
    </div>
  );
}
