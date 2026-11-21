import { useEffect, useMemo, useRef, useState } from "react";
import type { CorpusSample } from "../../api";

interface Props {
  samples: readonly CorpusSample[];
  onPick: (id: number) => void;
  onClose: () => void;
  notes?: ReadonlyMap<number, string>;
}

export function ComparePicker({
  samples,
  onPick,
  onClose,
  notes,
}: Props): JSX.Element {
  const [q, setQ] = useState("");
  const [cursorId, setCursorId] = useState<number | undefined>(undefined);
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
    const f = q.trim().toLowerCase();
    return f ? samples.filter((s) => s.name.toLowerCase().includes(f)) : samples;
  }, [samples, q]);
  useEffect(() => setCursorId(undefined), [q]);
  const cursor = Math.max(0, shown.findIndex((s) => s.id === cursorId));
  const moveTo = (i: number): void => setCursorId(shown[i]?.id);
  const pick = (s: CorpusSample | undefined): void => {
    if (!s) return;
    onPick(s.id);
    onClose();
  };
  return (
    <div
      ref={root}
      role="dialog"
      aria-label="Compare with a sample"
      data-keys-trap
      data-testid="workspace-compare-picker"
      className="absolute left-0 top-full z-20 mt-1 w-72 rounded-sm border border-hair bg-plate shadow-lg"
      onKeyDown={(e) => {
        if (e.key === "Escape") {
          e.preventDefault();
          onClose();
        } else if (e.key === "ArrowDown") {
          e.preventDefault();
          moveTo(Math.max(0, Math.min(shown.length - 1, cursor + 1)));
        } else if (e.key === "ArrowUp") {
          e.preventDefault();
          moveTo(Math.max(0, cursor - 1));
        } else if (e.key === "Enter") {
          e.preventDefault();
          pick(shown[cursor]);
        }
      }}
    >
      <input
        ref={input}
        type="search"
        value={q}
        onChange={(e) => setQ(e.target.value)}
        placeholder="Find a sample"
        aria-label="Find a sample to compare"
        className="w-full border-b border-hair bg-plate px-2 py-1 text-body placeholder:text-ink-faint focus:outline-none"
      />
      <ul className="max-h-64 overflow-y-auto py-1">
        {shown.map((s, i) => (
          <li key={s.id}>
            <button
              type="button"
              data-sample-id={s.id}
              aria-label={notes?.get(s.id) ? `${s.name}, ${notes.get(s.id)}` : undefined}
              aria-current={i === cursor ? "true" : undefined}
              onMouseEnter={() => setCursorId(s.id)}
              onClick={() => pick(s)}
              className={"flex w-full items-baseline gap-2 px-3 py-1 text-left text-body " + (i === cursor ? "bg-row-cursor" : "")}
            >
              <span className="min-w-0 flex-1 truncate">{s.name}</span>
              {notes?.get(s.id) && <span className="shrink-0 text-caption">{notes.get(s.id)}</span>}
            </button>
          </li>
        ))}
        {shown.length === 0 && <li className="px-3 py-2 text-caption">No samples match.</li>}
      </ul>
    </div>
  );
}
