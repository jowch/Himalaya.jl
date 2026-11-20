import { CUSTOM_PHASES, type CustomDraft } from "../indexing";
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

export function CustomIndexPanel({ draft, min, max, landed, total, disabled, onChange, onAdd, onCancel }: Props): JSX.Element {
  const param = draft.phase === "Lamellar" ? "d" : "a";
  return (
    <div
      className="flex flex-col gap-3"
      data-testid="workspace-custom-index"
      data-keys-trap
      onKeyDown={(e) => {
        if (e.key === "Escape") {
          e.preventDefault();
          onCancel();
        } else if (e.key === "Enter" && !disabled) {
          e.preventDefault();
          onAdd();
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
          value={draft.lattice}
          aria-label={`Lattice ${param} value`}
          onChange={(e) => {
            const v = Number(e.target.value);
            if (Number.isFinite(v) && v > 0) onChange({ ...draft, lattice: v });
          }}
          className="w-20 rounded-sm border border-hair bg-plate px-1 py-0.5 text-right text-body"
        />
        <span className="text-ink-soft">Å</span>
      </label>
      <p className="text-caption" data-testid="workspace-custom-landed">
        <span className="font-semibold text-ink">
          {landed} of {total}
        </span>{" "}
        reflections land on a peak
      </p>
      <p className="text-caption">
        <Kbd>←</Kbd> <Kbd>→</Kbd> nudge the lattice, Shift for larger steps
      </p>
      <div className="flex gap-2">
        <button
          type="button"
          onClick={onAdd}
          disabled={disabled}
          className="flex flex-1 items-center justify-between rounded-sm bg-accent px-3 py-1.5 text-body font-semibold text-paper disabled:opacity-50"
        >
          <span>
            Add {draft.phase} {draft.lattice.toFixed(1)} Å
          </span>
          <Kbd>Enter</Kbd>
        </button>
        <button type="button" onClick={onCancel} className="inline-flex items-center gap-1 rounded-sm border border-hair px-3 py-1.5 text-body">
          Cancel <Kbd>Esc</Kbd>
        </button>
      </div>
    </div>
  );
}
