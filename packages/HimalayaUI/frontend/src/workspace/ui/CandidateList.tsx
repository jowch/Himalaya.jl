import type { IndexEntry } from "../../api";
import { PhaseChip } from "../../print/ui/PhaseChip";
import type { CallRow } from "../bonnetGroups";
import { Kbd } from "./Kbd";
import { Note } from "./Note";

interface Props {
  rows: readonly CallRow[];
  assigned: ReadonlySet<number>;
  disabled: boolean;
  formFactor: boolean;
  onToggle: (ix: IndexEntry) => void;
  onToggleGroup: (members: readonly IndexEntry[]) => void;
  onFormFactor: () => void;
  onHover: (ids: readonly number[] | undefined) => void;
}

const ROW = "grid w-full items-center gap-2 rounded-sm border px-2 py-[5px] text-left text-base disabled:opacity-50 ";
const look = (on: boolean) => (on ? "border-ink-soft bg-paper-sunk" : "border-transparent hover:border-hair-strong");

export function CandidateList({ rows, assigned, disabled, formFactor, onToggle, onToggleGroup, onFormFactor, onHover }: Props): JSX.Element {
  return (
    <div className="flex flex-col gap-1">
      {rows.length === 0 && <Note>No candidate indices for this frame.</Note>}
      {rows.length > 0 && <Rows rows={rows} assigned={assigned} disabled={disabled} onToggle={onToggle} onToggleGroup={onToggleGroup} onHover={onHover} />}
      <button
        type="button"
        aria-pressed={formFactor}
        disabled={disabled}
        onClick={onFormFactor}
        data-testid="workspace-form-factor"
        className={ROW + "grid-cols-[auto_1fr] " + look(formFactor)}
      >
        <Kbd>0</Kbd>
        <span>No phase (form factor)</span>
      </button>
    </div>
  );
}

function Key({ n }: { n: number }): JSX.Element {
  return n <= 9 ? <Kbd>{n}</Kbd> : <span aria-hidden="true" className="invisible"><Kbd>9</Kbd></span>;
}

function Rows({ rows, assigned, disabled, onToggle, onToggleGroup, onHover }: Omit<Props, "formFactor" | "onFormFactor">): JSX.Element {
  let n = 0;
  const candidateRow = (c: IndexEntry, inGroup: boolean) => {
    const on = assigned.has(c.id);
    const key = ++n;
    return (
      <li key={c.id} className={inGroup ? "pl-4" : undefined}>
        <button
          type="button"
          aria-pressed={on}
          disabled={disabled}
          onClick={() => onToggle(c)}
          onMouseEnter={() => onHover([c.id])}
          onFocus={() => onHover([c.id])}
          className={ROW + "grid-cols-[auto_auto_1fr_auto] " + look(on)}
        >
          <Key n={key} />
          <PhaseChip phase={c.phase} variant="solid" />
          <span className="text-ink-soft">{c.lattice_d !== null ? `${c.lattice_d.toFixed(1)} Å` : ""}</span>
          <span className="text-right text-ink-faint">{c.kind === "speculative" ? "Custom" : c.score !== null ? c.score.toFixed(2) : ""}</span>
        </button>
      </li>
    );
  };
  return (
    <ol className="flex flex-col gap-1" data-testid="workspace-candidates" onMouseLeave={() => onHover(undefined)}>
      {rows.map((row) => {
        if (row.kind === "one") return candidateRow(row.candidate, false);
        const ids = row.members.map((m) => m.id);
        const count = ids.filter((id) => assigned.has(id)).length;
        const pressed = count === ids.length ? true : count === 0 ? false : "mixed";
        const key = ++n;
        return (
          <li key={`group-${ids.join("-")}`}>
            <button
              type="button"
              aria-pressed={pressed}
              disabled={disabled}
              onClick={() => onToggleGroup(row.members)}
              onMouseEnter={() => onHover(ids)}
              onFocus={() => onHover(ids)}
              data-testid="workspace-bonnet-group"
              className={ROW + "grid-cols-[auto_1fr_auto] " + (pressed === "mixed" ? "border-dashed border-ink-soft" : look(pressed === true))}
            >
              <Key n={key} />
              <span>Bonnet group</span>
              <span className="text-right text-ink-faint">{row.peaks} peaks</span>
            </button>
            <ol className="mt-1 flex flex-col gap-1">{row.members.map((m) => candidateRow(m, true))}</ol>
          </li>
        );
      })}
    </ol>
  );
}
