import { useEffect, useMemo, useRef, useState } from "react";
import type { CorpusSample } from "../../api";

interface Props {
  samples: readonly CorpusSample[];
  onPick: (id: number) => void;
  onClose: () => void;
}

export function ComparePicker({ samples, onPick, onClose }: Props): JSX.Element {
  const [q, setQ] = useState("");
  const [cursor, setCursor] = useState(0);
  const input = useRef<HTMLInputElement | null>(null);
  useEffect(() => input.current?.focus(), []);
  const shown = useMemo(() => {
    const f = q.trim().toLowerCase();
    return f ? samples.filter((s) => s.name.toLowerCase().includes(f)) : samples;
  }, [samples, q]);
  useEffect(() => setCursor(0), [q]);
  const pick = (s: CorpusSample | undefined): void => {
    if (!s) return;
    onPick(s.id);
    onClose();
  };
  return (
    <div
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
          setCursor((c) => Math.min(shown.length - 1, c + 1));
        } else if (e.key === "ArrowUp") {
          e.preventDefault();
          setCursor((c) => Math.max(0, c - 1));
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
              aria-current={i === cursor ? "true" : undefined}
              onMouseEnter={() => setCursor(i)}
              onClick={() => pick(s)}
              className={"w-full px-3 py-1 text-left text-body " + (i === cursor ? "bg-row-cursor" : "")}
            >
              {s.name}
            </button>
          </li>
        ))}
        {shown.length === 0 && <li className="px-3 py-2 text-caption">No samples match.</li>}
      </ul>
    </div>
  );
}
