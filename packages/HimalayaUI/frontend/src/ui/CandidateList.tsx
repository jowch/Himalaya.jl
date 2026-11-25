import { useId, useRef, useState } from "react";
import type { IndexEntry } from "../api";
import { PhaseChip } from "./PhaseChip";
import type { CallRow } from "../workspace/bonnetGroups";
import { formatCurvature, gaussianCurvature } from "../workspace/curvature";
import { Kbd } from "./Kbd";
import { Note } from "./Note";
import { Tooltip } from "./Tooltip";
import { CURVATURE, NOT_FOUND, SCORE } from "./explanations";
import { Explanation, Term } from "./Term";

interface Props {
  rows: readonly CallRow[];
  assigned: ReadonlySet<number>;
  unsupported?: ReadonlySet<number>;
  disabled: boolean;
  formFactor: boolean;
  onToggle: (ix: IndexEntry) => void;
  onToggleGroup: (members: readonly IndexEntry[]) => void;
  onFormFactor: () => void;
  onHover: (ids: readonly number[] | undefined) => void;
}

const NONE: ReadonlySet<number> = new Set();

const ROW = "grid w-full items-center rounded-sm border px-2 text-left text-base disabled:opacity-50 ";
const ROW_LOOK = "border-transparent hover:bg-ink/6";

export function CandidateList({ rows, assigned, unsupported = NONE, disabled, formFactor, onToggle, onToggleGroup, onFormFactor, onHover }: Props): JSX.Element {
  const [pointerFocus, setPointerFocus] = useState(false);
  const pointerDown = useRef(false);
  return (
    <div
      className="flex flex-col gap-1"
      role="group"
      aria-label="Phase choices"
      data-call-choices=""
      data-pointer-focus={pointerFocus ? "" : undefined}
      onPointerDown={() => {
        pointerDown.current = true;
        setPointerFocus(true);
      }}
      onPointerUp={() => {
        pointerDown.current = false;
      }}
      onPointerCancel={() => {
        pointerDown.current = false;
      }}
      onFocus={() => {
        setPointerFocus(pointerDown.current);
        pointerDown.current = false;
      }}
      onBlur={(e) => {
        if (e.currentTarget.contains(e.relatedTarget)) return;
        onHover(undefined);
        setPointerFocus(false);
      }}
    >
      {rows.every((r) => r.kind === "one" && unsupported.has(r.candidate.id)) && <Note>No candidate indices for this image.</Note>}
      {rows.length > 0 && <Rows rows={rows} assigned={assigned} unsupported={unsupported} disabled={disabled} onToggle={onToggle} onToggleGroup={onToggleGroup} onHover={onHover} />}
      <button
        type="button"
        aria-pressed={formFactor}
        disabled={disabled}
        onClick={onFormFactor}
        onFocus={() => onHover(undefined)}
        data-testid="workspace-form-factor"
        className={ROW + "grid-cols-[auto_1fr] gap-2 py-[5px] " + ROW_LOOK}
      >
        <Kbd pressed={formFactor}>0</Kbd>
        <span>No phase (form factor)</span>
      </button>
    </div>
  );
}

function Key({ n, pressed }: { n: number; pressed: boolean | "mixed" }): JSX.Element {
  return n <= 9 ? <Kbd pressed={pressed}>{n}</Kbd> : <span aria-hidden="true" data-blank-key=""><Kbd pressed={pressed}>{"\u00a0"}</Kbd></span>;
}

function Rows({ rows, assigned, unsupported, disabled, onToggle, onToggleGroup, onHover }: Omit<Props, "formFactor" | "onFormFactor" | "unsupported"> & { unsupported: ReadonlySet<number> }): JSX.Element {
  const base = useId();
  let n = 0;
  const candidateRow = (c: IndexEntry, inGroup: boolean) => {
    const on = assigned.has(c.id);
    const missing = unsupported.has(c.id);
    const key = ++n;
    const lattice = c.lattice_d !== null ? `${c.lattice_d.toFixed(1)} Å` : "";
    const score = missing || c.kind === "speculative" ? null : c.score;
    const why = missing ? NOT_FOUND : score !== null ? SCORE : undefined;
    const whyId = `${base}-${c.id}`;
    return (
      <li key={c.id} className={inGroup ? "pl-2" : undefined}>
        <button
          type="button"
          aria-pressed={on}
          disabled={disabled}
          onClick={() => onToggle(c)}
          onMouseEnter={() => onHover([c.id])}
          onMouseMove={() => onHover([c.id])}
          onFocus={() => onHover([c.id])}
          data-unsupported={missing ? "" : undefined}
          {...(why ? { "aria-describedby": whyId } : {})}
          className={ROW + "grid-cols-[auto_auto_1fr_auto] py-[5px] " + (inGroup ? "gap-1 " : "gap-1.5 ") + ROW_LOOK}
        >
          <Key n={key} pressed={on} />
          <PhaseChip phase={c.phase} variant="text" size="md" />
          <Tooltip content={lattice} overflowOnly>
            <span className="min-w-0 truncate whitespace-nowrap text-ink-soft">{lattice}</span>
          </Tooltip>
          <span className="text-right text-ink-faint">
            {missing ? (
              <Term explain={NOT_FOUND} nested>Not found</Term>
            ) : score !== null ? (
              <Term explain={SCORE} nested>{score.toFixed(2)}</Term>
            ) : c.kind === "speculative" ? "Custom" : ""}
          </span>
          {why && <Explanation id={whyId} texts={[why]} />}
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
        const kappas = row.members.map((m) => gaussianCurvature(m.phase, m.lattice_d)).filter((k): k is number => k !== null);
        const shared = kappas.length > 0 ? kappas.reduce((x, y) => x + y, 0) / kappas.length : null;
        return (
          <li key={`group-${ids.join("-")}`}>
            <button
              type="button"
              aria-pressed={pressed}
              disabled={disabled}
              onClick={() => onToggleGroup(row.members)}
              onMouseEnter={() => onHover(ids)}
              onMouseMove={() => onHover(ids)}
              onFocus={() => onHover(ids)}
              data-testid="workspace-bonnet-group"
              {...(shared !== null ? { "aria-describedby": `${base}-group-${ids.join("-")}` } : {})}
              className={ROW + "grid-cols-[auto_1fr_auto] gap-1.5 py-[5px] " + ROW_LOOK}
            >
              <Key n={key} pressed={pressed} />
              <span className="min-w-0 truncate whitespace-nowrap">
                Bonnet group
                {shared !== null && <Term explain={CURVATURE} nested className="ml-1.5 text-xs text-ink-soft"><span data-testid="workspace-bonnet-curvature">{formatCurvature(shared)}</span></Term>}
              </span>
              <span className="text-right text-ink-faint">{row.peaks} peaks</span>
              {shared !== null && <Explanation id={`${base}-group-${ids.join("-")}`} texts={[CURVATURE]} />}
            </button>
            <ol className="mt-1 flex flex-col gap-1">{row.members.map((m) => candidateRow(m, true))}</ol>
          </li>
        );
      })}
    </ol>
  );
}
