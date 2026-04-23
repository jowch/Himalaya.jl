import { useState, useMemo } from "react";
import type { Sample } from "../api";

export interface SampleListProps {
  samples: Sample[];
  activeId: number | undefined;
  onSelect: (id: number) => void;
}

type SampleStatus = "unanalyzed" | "candidates" | "confirmed";

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
    <div className="sample-list-wrap">
      <div className="sample-list-header">
        <input
          className="sample-filter"
          type="text"
          placeholder="Filter samples…"
          value={filter}
          onChange={(e) => setFilter(e.target.value)}
        />
      </div>
      <ul className="sample-list">
        {filtered.map((s) => {
          const status: SampleStatus = "unanalyzed";
          return (
            <li
              key={s.id}
              className={"sample-row" + (s.id === activeId ? " active" : "")}
              data-sample-id={s.id}
              data-status={status}
              onClick={() => onSelect(s.id)}
            >
              <span className={`status-dot status-${status}`} />
              <span className="sample-label">{s.label ?? ""}</span>
              <span className="sample-name muted">{s.name ?? ""}</span>
            </li>
          );
        })}
      </ul>
    </div>
  );
}
