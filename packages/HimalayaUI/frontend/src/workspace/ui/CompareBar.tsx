import { Kbd } from "./Kbd";
import type { CompareScale } from "../compare";
import { COMPARE_DASH } from "./WorkspacePlot";

export interface CompareChip {
  id: number;
  name: string;
  color: string;
  scale: CompareScale | undefined;
  bridged?: boolean;
}

interface Props {
  chips: readonly CompareChip[];
  subtractId: number | undefined;
  hidden: boolean;
  onSubtract: (id: number) => void;
  onRemove: (id: number) => void;
  onAdd: () => void;
}

export function CompareBar({ chips, subtractId, hidden, onSubtract, onRemove, onAdd }: Props): JSX.Element {
  return (
    <div className="flex min-w-0 flex-nowrap items-center gap-2" data-testid="workspace-compare-bar">
      {chips.map((c) => {
        const sub = c.id === subtractId;
        return (
          <span
            key={c.id}
            data-compare-chip={c.id}
            style={{ color: c.color }}
            className={
              "inline-flex min-w-0 items-center gap-1.5 rounded-sm border bg-plate py-0.5 pl-2 pr-1 text-caption text-ink " +
              (sub ? "border-ink-soft" : "border-hair") +
              (hidden ? " opacity-45" : "")
            }
          >
            <svg width={18} height={4} aria-hidden="true" className="shrink-0">
              <line x1={0} x2={18} y1={2} y2={2} stroke="currentColor" strokeWidth={2} strokeDasharray={COMPARE_DASH} />
            </svg>
            <button
              type="button"
              aria-pressed={sub}
              {...(hidden ? { "aria-label": `${c.name}, hidden` } : {})}
              title={sub ? "Stop subtracting" : "Subtract from the current trace"}
              onClick={() => onSubtract(c.id)}
              className="min-w-0 truncate text-ink"
            >
              {c.name}
            </button>
            {c.scale !== undefined && (
              <span
                className="shrink-0 whitespace-nowrap text-ink-faint"
                data-scale-source={c.scale.source}
                title={c.scale.source === "monitor" ? "Scaled by the transmitted beam" : "No beam monitors, so matched at high q"}
              >
                ×{c.scale.scale.toPrecision(3)}
                {c.scale.source === "matched" && " matched"}
              </span>
            )}
            {c.bridged && (
              <span className="shrink-0 text-ink-faint" data-bridged-chip>
                bridged
              </span>
            )}
            <button
              type="button"
              aria-label={c.bridged ? `Remove bridge to ${c.name}` : `Remove ${c.name}`}
              title={c.bridged ? "Take this sample out of the bridge" : undefined}
              onClick={() => onRemove(c.id)}
              className="shrink-0 px-1 text-ink-faint hover:text-ink"
            >
              ×
            </button>
          </span>
        );
      })}
      <button type="button" onClick={onAdd} className="inline-flex shrink-0 items-center gap-1 whitespace-nowrap rounded-sm border border-dashed border-hair-strong px-2 py-0.5 text-caption text-ink-soft hover:text-ink">
        + Compare <Kbd>C</Kbd>
      </button>
    </div>
  );
}
