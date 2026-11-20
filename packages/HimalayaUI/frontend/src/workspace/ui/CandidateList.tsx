import type { IndexEntry } from "../../api";
import { PhaseChip } from "../../print/ui/PhaseChip";
import { Kbd } from "./Kbd";
import { Note } from "./Note";

interface Props {
  candidates: readonly IndexEntry[];
  assigned: ReadonlySet<number>;
  disabled: boolean;
  onToggle: (ix: IndexEntry) => void;
  onHover: (id: number | undefined) => void;
}

export function CandidateList({ candidates, assigned, disabled, onToggle, onHover }: Props): JSX.Element {
  if (candidates.length === 0) return <Note>No candidate indices for this exposure.</Note>;
  return (
    <ol className="flex flex-col gap-0.5" data-testid="workspace-candidates" onMouseLeave={() => onHover(undefined)}>
      {candidates.map((c, i) => {
        const on = assigned.has(c.id);
        return (
          <li key={c.id}>
            <button
              type="button"
              aria-pressed={on}
              disabled={disabled}
              onClick={() => onToggle(c)}
              onMouseEnter={() => onHover(c.id)}
              onFocus={() => onHover(c.id)}
              className={
                "grid w-full grid-cols-[auto_auto_1fr_auto] items-center gap-2 rounded-sm border px-2 py-1 text-left text-body disabled:opacity-50 " +
                (on ? "border-accent bg-accent-wash" : "border-transparent hover:bg-paper-sunk")
              }
            >
              <Kbd>{i + 1}</Kbd>
              <PhaseChip phase={c.phase} variant="solid" />
              <span>{c.lattice_d !== null ? `${c.lattice_d.toFixed(1)} Å` : ""}</span>
              <span className="text-ink-soft">{c.kind === "speculative" ? "Custom" : c.score !== null ? c.score.toFixed(2) : ""}</span>
            </button>
          </li>
        );
      })}
    </ol>
  );
}
