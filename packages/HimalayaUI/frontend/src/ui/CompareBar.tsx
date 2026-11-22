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
  setup?: string | undefined;
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
  full?: boolean;
  trailing?: ReactNode;
}

export function CompareBar({ chips, subtractId, hidden, modKey, onToggleShown, onSubtract, onRemove, onAdd, full = false, trailing }: Props): JSX.Element {
  const perRow = Math.max(1, Math.ceil((chips.length + (full ? 0 : 1) + (trailing ? 1 : 0)) / 2));
  const basis = `calc((100% - ${perRow - 1} * 0.5rem) / ${perRow})`;
  return (
    <div className="flex min-w-0 items-start gap-2" data-testid="workspace-compare-bar">
      <div data-compare-chips className="flex h-[calc(2*var(--compare-row)+0.25rem)] min-w-0 flex-1 flex-wrap content-start items-center gap-x-2 gap-y-1 [--compare-row:1.5rem]">
        {chips.map((c) => {
          const sub = c.id === subtractId;
          const off = hidden || c.hidden === true;
          const dim = off ? " opacity-45" : "";
          return (
            <span
              key={c.id}
              data-compare-chip={c.id}
              style={{ color: c.color, flexBasis: basis }}
              className={
                "inline-flex h-(--compare-row) min-w-min max-w-max items-center gap-1.5 rounded-sm border bg-plate py-0.5 pl-2 pr-1 text-caption text-ink " +
                (sub ? "border-accent" : "border-hair")
              }
            >
              <svg width={18} height={4} aria-hidden="true" className={"shrink-0" + dim}>
                <line x1={0} x2={18} y1={2} y2={2} stroke="currentColor" strokeWidth={2} strokeDasharray={COMPARE_DASH} />
              </svg>
              <button
                type="button"
                aria-pressed={!off}
                {...(off || c.setup ? { "aria-label": [c.name, c.setup, off ? "hidden" : undefined].filter(Boolean).join(", ") } : {})}
                aria-description={sub ? "Subtracted from the current trace" : undefined}
                data-subtracted={sub ? "" : undefined}
                title={`${c.setup ? `From ${c.setup}. ` : ""}Click to ${off ? "show" : "hide"}. ${modKey}-click to ${sub ? "stop subtracting" : "subtract from the current trace"}.`}
                onClick={(e) => {
                  if (e.metaKey || e.ctrlKey) onSubtract(c.id);
                  else onToggleShown(c.id);
                }}
                className={"grid min-w-min grow grid-cols-[minmax(4ch,1fr)_auto] items-center gap-1 text-left text-ink" + dim}
              >
                <span className="min-w-0 truncate">{c.name}</span>
                {c.setup && (
                  <svg data-chip-setup width={10} height={10} viewBox="0 0 10 10" aria-hidden="true" className="text-ink-faint">
                    <rect x={0.5} y={3.5} width={6} height={6} rx={1} fill="none" stroke="currentColor" />
                    <path d="M3.5 1.5 a1 1 0 0 1 1 -1 h4 a1 1 0 0 1 1 1 v4 a1 1 0 0 1 -1 1" fill="none" stroke="currentColor" />
                  </svg>
                )}
              </button>
              {c.scale !== undefined && (
                <span
                  className={"grid min-w-min grid-cols-[auto_minmax(0,1fr)] items-baseline whitespace-nowrap text-ink-faint" + dim}
                  data-scale-source={c.scale.source}
                  title={c.scale.source === "monitor" ? "Scaled by the transmitted beam" : "No beam monitors, so matched at high q"}
                >
                  <span data-scale-value className="shrink-0">×{c.scale.scale.toPrecision(3)}</span>
                  {c.scale.source === "matched" && <span className="min-w-0 truncate">{"\u00a0matched"}</span>}
                </span>
              )}
              <button
                type="button"
                aria-label={`Remove ${[c.name, c.setup].filter(Boolean).join(", ")}`}
                onClick={() => onRemove(c.id)}
                className="shrink-0 px-1 text-ink-faint hover:text-ink"
              >
                ×
              </button>
            </span>
          );
        })}
        {!full && (
          <button type="button" onClick={onAdd} className="inline-flex h-(--compare-row) shrink-0 items-center gap-1 whitespace-nowrap rounded-sm border border-dashed border-hair-strong px-2 py-0.5 text-caption text-ink-soft hover:text-ink">
            + Compare <Kbd>C</Kbd>
          </button>
        )}
        {trailing && <div className="ml-auto flex h-(--compare-row) shrink-0 items-center">{trailing}</div>}
      </div>
    </div>
  );
}
