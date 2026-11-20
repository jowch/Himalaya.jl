import { Kbd } from "./Kbd";

export interface CompareChip {
  id: number;
  name: string;
  dash: string;
}

interface Props {
  chips: readonly CompareChip[];
  subtractId: number | undefined;
  hidden: boolean;
  logY: boolean;
  onToggleLog: () => void;
  onSubtract: (id: number) => void;
  onRemove: (id: number) => void;
  onAdd: () => void;
}

export function CompareBar({ chips, subtractId, hidden, logY, onToggleLog, onSubtract, onRemove, onAdd }: Props): JSX.Element {
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
            {hidden ? "Comparisons hidden" : "Comparisons shown"} <Kbd>R</Kbd>
          </span>
        )}
        <button type="button" onClick={onToggleLog} aria-pressed={logY} className="inline-flex items-center gap-1 hover:text-ink">
          {logY ? "Log" : "Linear"} <Kbd>L</Kbd>
        </button>
      </span>
    </div>
  );
}
