import { useState, useMemo } from "react";
import type { Sample } from "../api";

export interface SampleListProps {
  samples: Sample[];
  activeId: number | undefined;
  onSelect: (id: number) => void;
}

type SampleStatus = "unanalyzed" | "candidates" | "confirmed";

const statusDot: Record<SampleStatus, string> = {
  unanalyzed: "bg-fg-dim",
  candidates: "bg-warning",
  confirmed:  "bg-success",
};

function matches(s: Sample, filter: string): boolean {
  if (!filter) return true;
  const n = filter.toLowerCase();
  if (s.label?.toLowerCase().includes(n)) return true;
  if (s.name?.toLowerCase().includes(n))  return true;
  return s.tags.some((t) =>
    t.key.toLowerCase().includes(n) || t.value.toLowerCase().includes(n),
  );
}

export function SampleList({ samples, activeId, onSelect }: SampleListProps): JSX.Element {
  const [filter, setFilter] = useState("");
  const filtered = useMemo(() => samples.filter((s) => matches(s, filter)), [samples, filter]);

  return (
    <div className="flex flex-col h-full">
      <div className="p-2 border-b border-border">
        <input
          className="w-full bg-bg-elevated border border-border rounded-md px-2 py-1 focus:outline focus:outline-1 focus:outline-accent focus:border-accent"
          type="text"
          placeholder="Filter samples…"
          value={filter}
          onChange={(e) => setFilter(e.target.value)}
        />
      </div>
      <ul className="list-none overflow-y-auto flex-1">
        {filtered.map((s) => {
          const status: SampleStatus = "unanalyzed";
          const active = s.id === activeId;
          return (
            <li
              key={s.id}
              className={
                "grid grid-cols-[12px_auto_1fr] items-center gap-2 py-1.5 px-4 cursor-pointer border-l-2 " +
                (active ? "bg-bg-elevated border-accent" : "border-transparent hover:bg-bg-hover")
              }
              data-sample-id={s.id}
              data-active={active}
              data-status={status}
              onClick={() => onSelect(s.id)}
            >
              <span className={`w-2 h-2 rounded-full ${statusDot[status]}`} />
              <span className="font-medium" data-testid="sample-label">{s.label ?? ""}</span>
              <span className="text-fg-muted">{s.name ?? ""}</span>
            </li>
          );
        })}
      </ul>
    </div>
  );
}
