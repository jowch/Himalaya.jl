import type { ReactNode } from "react";
import { Kbd } from "./Kbd";
import type { CompareScale } from "../workspace/compare";
import { COMPARE_DASH } from "./WorkspacePlot";

export interface CompareChip {
  id: number;
  name: string;
  color: string;
  scale: CompareScale | undefined;
  hidden?: boolean;
}

interface Props {
  chips: readonly CompareChip[];
  subtractId: number | undefined;
  hidden: boolean;
  modKey: string;
  onToggleShown: (id: number) => void;
  onSubtract: (id: number) => void;
  onRemove: (id: number) => void;
  onAdd: () => void;
  trailing?: ReactNode;
}

export function CompareBar({ chips, subtractId, hidden, modKey, onToggleShown, onSubtract, onRemove, onAdd, trailing }: Props): JSX.Element {
  return (
    <div className="flex min-w-0 flex-nowrap items-center gap-2" data-testid="workspace-compare-bar">
      {chips.map((c) => {
        const sub = c.id === subtractId;
        const off = hidden || c.hidden === true;
        const dim = off ? " opacity-45" : "";
        return (
          <span
            key={c.id}
            data-compare-chip={c.id}
            style={{ color: c.color }}
            className={
              "inline-flex min-w-0 items-center gap-1.5 rounded-sm border bg-plate py-0.5 pl-2 pr-1 text-caption text-ink " +
              (sub ? "border-accent" : "border-hair")
            }
          >
            <svg width={18} height={4} aria-hidden="true" className={"shrink-0" + dim}>
              <line x1={0} x2={18} y1={2} y2={2} stroke="currentColor" strokeWidth={2} strokeDasharray={COMPARE_DASH} />
            </svg>
            <button
              type="button"
              aria-pressed={!off}
              {...(off ? { "aria-label": `${c.name}, hidden` } : {})}
              aria-description={sub ? "Subtracted from the current trace" : undefined}
              data-subtracted={sub ? "" : undefined}
              title={`Click to ${off ? "show" : "hide"}. ${modKey}-click to ${sub ? "stop subtracting" : "subtract from the current trace"}.`}
              onClick={(e) => {
                if (e.metaKey || e.ctrlKey) onSubtract(c.id);
                else onToggleShown(c.id);
              }}
              className={"min-w-0 truncate text-ink" + dim}
            >
              {c.name}
            </button>
            {c.scale !== undefined && (
              <span
                className={"shrink-0 whitespace-nowrap text-ink-faint" + dim}
                data-scale-source={c.scale.source}
                title={c.scale.source === "monitor" ? "Scaled by the transmitted beam" : "No beam monitors, so matched at high q"}
              >
                ×{c.scale.scale.toPrecision(3)}
                {c.scale.source === "matched" && " matched"}
              </span>
            )}
            <button
              type="button"
              aria-label={`Remove ${c.name}`}
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
      {trailing && <div className="ml-auto shrink-0">{trailing}</div>}
    </div>
  );
}
