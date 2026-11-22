import { useEffect, useRef, useState } from "react";
import { CUSTOM_PHASES, type CustomDraft } from "../workspace/indexing";
import { Kbd } from "./Kbd";

interface Props {
  draft: CustomDraft;
  min: number;
  max: number;
  landed: number;
  total: number;
  disabled: boolean;
  onChange: (d: CustomDraft) => void;
  onAdd: () => void;
  onCancel: () => void;
}

export const NUDGE_STEP = 0.5;
export const NUDGE_STEP_LARGE = 5;

function parseLattice(text: string, min: number, max: number): number | null {
  if (text.trim() === "") return null;
  const v = Number(text);
  return Number.isFinite(v) && v >= min && v <= max ? v : null;
}

export function CustomIndexPanel({ draft, min, max, landed, total, disabled, onChange, onAdd, onCancel }: Props): JSX.Element {
  const param = draft.phase === "Lamellar" ? "d" : "a";
  const [text, setText] = useState(draft.lattice.toFixed(1));
  const editing = useRef(false);
  useEffect(() => {
    if (!editing.current) setText(draft.lattice.toFixed(1));
  }, [draft.lattice]);
  const typed = parseLattice(text, min, max);
  const invalid = typed === null;
  const canAdd = !disabled && !invalid;

  const settle = (): void => {
    editing.current = false;
    const v = Number(text);
    const next = Number.isFinite(v) && text.trim() !== "" ? Math.min(max, Math.max(min, v)) : draft.lattice;
    onChange({ ...draft, lattice: next });
    setText(next.toFixed(1));
  };

  return (
    <div
      className="flex flex-col gap-3"
      data-testid="workspace-custom-index"
      data-keys-trap
      onKeyDown={(e) => {
        const target = e.target as HTMLElement;
        const inField = target.tagName === "INPUT";
        if (e.key === "Escape") {
          e.preventDefault();
          onCancel();
        } else if (e.key === "Enter" && target.tagName !== "BUTTON") {
          e.preventDefault();
          if (invalid) settle();
          else if (canAdd) onAdd();
        } else if ((e.key === "ArrowLeft" || e.key === "ArrowRight") && !inField) {
          e.preventDefault();
          const step = (e.key === "ArrowRight" ? 1 : -1) * (e.shiftKey ? NUDGE_STEP_LARGE : NUDGE_STEP);
          onChange({ ...draft, lattice: draft.lattice + step });
        }
      }}
    >
      <div className="flex flex-wrap gap-1" role="group" aria-label="Symmetry">
        {CUSTOM_PHASES.map((p) => (
          <button
            key={p}
            type="button"
            aria-pressed={p === draft.phase}
            onClick={() => onChange({ ...draft, phase: p })}
            className={
              "rounded-sm border px-2 py-0.5 text-caption " +
              (p === draft.phase ? "border-ink bg-ink text-paper" : "border-hair text-ink-soft hover:border-hair-strong")
            }
          >
            {p}
          </button>
        ))}
      </div>
      <label className="flex items-center gap-2 text-body">
        <span className="text-ink-soft">{param}</span>
        <input
          type="range"
          min={min}
          max={max}
          step={0.1}
          value={draft.lattice}
          aria-label={`Lattice ${param} in ångström`}
          onChange={(e) => onChange({ ...draft, lattice: Number(e.target.value) })}
          className="min-w-0 flex-1 accent-accent"
        />
        <input
          type="number"
          min={min}
          max={max}
          step={0.1}
          value={text}
          aria-label={`Lattice ${param} value`}
          aria-invalid={invalid}
          onFocus={() => {
            editing.current = true;
          }}
          onChange={(e) => {
            editing.current = true;
            setText(e.target.value);
            const v = parseLattice(e.target.value, min, max);
            if (v !== null) onChange({ ...draft, lattice: v });
          }}
          onBlur={settle}
          className="w-20 rounded-sm border border-hair bg-plate px-1 py-0.5 text-right text-body aria-[invalid=true]:border-error"
        />
        <span className="text-ink-soft">Å</span>
      </label>
      {invalid ? (
        <p className="text-caption text-error" data-testid="workspace-custom-invalid">
          Enter a value from {min.toFixed(1)} to {max.toFixed(1)} Å.
        </p>
      ) : total === 0 ? (
        <p className="text-caption" data-testid="workspace-custom-landed">
          No reflections fall inside this trace&apos;s q range.
        </p>
      ) : (
        <p className="text-caption" data-testid="workspace-custom-landed">
          <span className="font-semibold text-ink">
            {landed} of {total}
          </span>{" "}
          reflections land on a peak
        </p>
      )}
      <div className="flex flex-col gap-1.5">
        <button
          type="button"
          onClick={onAdd}
          disabled={!canAdd}
          className="flex w-full items-center justify-between gap-2 rounded-sm bg-accent px-3 py-2 text-body font-semibold text-paper disabled:opacity-50"
        >
          <span className="min-w-0 truncate">
            Add {draft.phase} {(typed ?? draft.lattice).toFixed(1)} Å
          </span>
          <Kbd>Enter</Kbd>
        </button>
        <button type="button" onClick={onCancel} className="flex w-full items-center justify-between gap-2 rounded-sm border border-hair bg-paper px-3 py-1.5 text-body">
          <span>Cancel</span>
          <Kbd>Esc</Kbd>
        </button>
      </div>
    </div>
  );
}
