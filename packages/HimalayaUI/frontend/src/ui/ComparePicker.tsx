import { useEffect, useMemo, useRef, useState } from "react";
import type { CorpusSample } from "../api";
import { useListHighlight } from "./useListHighlight";

interface Props {
  samples: readonly CorpusSample[];
  onPick: (id: number) => void;
  onClose: () => void;
  notes?: ReadonlyMap<number, string>;
  setups?: ReadonlyMap<number, string>;
  full?: string | undefined;
  title?: string;
  findLabel?: string;
  fill?: boolean;
}

export function ComparePicker({
  samples,
  onPick,
  onClose,
  notes,
  setups,
  full,
  title = "Compare with a sample",
  findLabel = "Find a sample to compare",
  fill = false,
}: Props): JSX.Element {
  const [q, setQ] = useState("");
  const input = useRef<HTMLInputElement | null>(null);
  const root = useRef<HTMLDivElement | null>(null);
  useEffect(() => input.current?.focus(), []);
  useEffect(() => {
    const onDown = (ev: MouseEvent): void => {
      if (root.current && ev.target instanceof Node && !root.current.contains(ev.target)) onClose();
    };
    document.addEventListener("mousedown", onDown, true);
    return () => document.removeEventListener("mousedown", onDown, true);
  }, [onClose]);
  const shown = useMemo(() => {
    if (full) return [];
    const f = q.trim().toLowerCase();
    return f ? samples.filter((s) => `${s.name} ${setups?.get(s.id) ?? ""}`.toLowerCase().includes(f)) : samples;
  }, [samples, q, full, setups]);
  const pick = (s: CorpusSample | undefined): void => {
    if (!s) return;
    onPick(s.id);
    onClose();
  };
  const shownIds = useMemo(() => shown.map((s) => s.id), [shown]);
  const highlight = useListHighlight({ count: shown.length, onCommit: (i) => pick(shown[i]), onEscape: onClose, keys: shownIds, resetKey: q });
  const cursor = highlight.at;
  return (
    <div
      ref={root}
      role="dialog"
      aria-label={title}
      data-keys-trap
      data-testid="workspace-compare-picker"
      className={"absolute top-full z-20 mt-1 rounded-sm border border-hair bg-plate shadow-lg " + (fill ? "inset-x-3.5" : "left-0 w-72")}
      onKeyDown={highlight.onKeyDown}
    >
      <input
        ref={input}
        type="search"
        value={q}
        onChange={(e) => setQ(e.target.value)}
        placeholder="Find a sample"
        aria-label={findLabel}
        className="w-full border-b border-hair bg-plate px-2 py-1 text-body placeholder:text-ink-faint focus:outline-none"
      />
      <ul className="max-h-64 overflow-y-auto py-1">
        {shown.map((s, i) => (
          <li key={s.id}>
            <button
              type="button"
              data-sample-id={s.id}
              aria-label={notes?.get(s.id) || setups?.get(s.id) ? [s.name, setups?.get(s.id), notes?.get(s.id)].filter(Boolean).join(", ") : undefined}
              tabIndex={-1}
              aria-current={i === cursor ? "true" : undefined}
              onMouseEnter={() => highlight.setAt(i)}
              onMouseDown={(e) => e.preventDefault()}
              onClick={() => pick(s)}
              className={"flex w-full items-baseline gap-2 px-3 py-1 text-left text-body " + (i === cursor ? "bg-row-cursor" : "")}
            >
              <span className="min-w-0 flex-1 truncate">{s.name}</span>
              {setups?.get(s.id) && <span className="max-w-24 shrink-0 truncate text-caption text-ink-faint">{setups.get(s.id)}</span>}
              {notes?.get(s.id) && <span className="shrink-0 text-caption">{notes.get(s.id)}</span>}
            </button>
          </li>
        ))}
        {shown.length === 0 && <li className="px-3 py-2 text-caption">{full ?? "No samples match."}</li>}
      </ul>
    </div>
  );
}
