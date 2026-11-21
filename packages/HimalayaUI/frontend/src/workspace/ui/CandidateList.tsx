import type { IndexEntry } from "../../api";
import { PhaseChip } from "../../print/ui/PhaseChip";
import { Kbd } from "./Kbd";
import { Note } from "./Note";

interface Props {
  candidates: readonly IndexEntry[];
  assigned: ReadonlySet<number>;
  disabled: boolean;
  formFactor: boolean;
  onToggle: (ix: IndexEntry) => void;
  onFormFactor: () => void;
  onHover: (id: number | undefined) => void;
}

export function CandidateList({ candidates, assigned, disabled, formFactor, onToggle, onFormFactor, onHover }: Props): JSX.Element {
  return (
    <div className="flex flex-col gap-1">
      {candidates.length === 0 && <Note>No candidate indices for this exposure.</Note>}
      {candidates.length > 0 && <Rows candidates={candidates} assigned={assigned} disabled={disabled} onToggle={onToggle} onHover={onHover} />}
      <button
        type="button"
        aria-pressed={formFactor}
        disabled={disabled}
        onClick={onFormFactor}
        data-testid="workspace-form-factor"
        className={
          "grid w-full grid-cols-[auto_1fr] items-center gap-2 rounded-sm border px-2 py-[5px] text-left text-base disabled:opacity-50 " +
          (formFactor ? "border-accent bg-accent-wash" : "border-transparent hover:border-hair-strong")
        }
      >
        <Kbd>0</Kbd>
        <span>No phase (form factor)</span>
      </button>
    </div>
  );
}

function Rows({ candidates, assigned, disabled, onToggle, onHover }: Omit<Props, "formFactor" | "onFormFactor">): JSX.Element {
  return (
    <ol className="flex flex-col gap-1" data-testid="workspace-candidates" onMouseLeave={() => onHover(undefined)}>
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
                "grid w-full grid-cols-[auto_auto_1fr_auto] items-center gap-2 rounded-sm border px-2 py-[5px] text-left text-base disabled:opacity-50 " +
                (on ? "border-accent bg-accent-wash" : "border-transparent hover:border-hair-strong")
              }
            >
              <Kbd>{i + 1}</Kbd>
              <PhaseChip phase={c.phase} variant="solid" />
              <span className="text-ink-soft">{c.lattice_d !== null ? `${c.lattice_d.toFixed(1)} Å` : ""}</span>
              <span className="text-right text-ink-faint">{c.kind === "speculative" ? "Custom" : c.score !== null ? c.score.toFixed(2) : ""}</span>
            </button>
          </li>
        );
      })}
    </ol>
  );
}
