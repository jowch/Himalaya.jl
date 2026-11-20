import { useEffect, useRef } from "react";
import type { CorpusSample } from "../../api";
import { sampleStatus } from "../model";

interface Props {
  samples: readonly CorpusSample[];
  currentId: number | undefined;
  filter: string;
  onFilter: (v: string) => void;
  onPick: (id: number) => void;
}

function StatusMark({ s }: { s: CorpusSample }): JSX.Element {
  const st = sampleStatus(s);
  if (st.kind === "phase") return <span data-testid="workspace-status" className="text-xs text-ink-soft">{st.phase}</span>;
  if (st.kind === "form_factor") return <span data-testid="workspace-status" className="text-xs text-ink-soft">Form factor</span>;
  return <span data-testid="workspace-status" role="img" aria-label="Unreviewed" title="Unreviewed" className="text-xs text-ink-faint">○</span>;
}

export function SampleList({ samples, currentId, filter, onFilter, onPick }: Props): JSX.Element {
  const currentRef = useRef<HTMLButtonElement | null>(null);
  useEffect(() => {
    currentRef.current?.scrollIntoView({ block: "nearest" });
  }, [currentId]);
  return (
    <>
      <div className="border-b border-hair p-2">
        <input
          id="workspace-filter"
          type="search"
          value={filter}
          onChange={(e) => onFilter(e.target.value)}
          placeholder="Filter  /"
          aria-label="Filter samples"
          className="w-full rounded-sm border border-hair bg-plate px-2 py-1 text-body placeholder:text-ink-faint focus:outline-2 focus:outline-accent"
        />
      </div>
      <ul className="min-h-0 flex-1 overflow-y-auto py-1" data-testid="workspace-samples">
        {samples.map((s) => {
          const cur = s.id === currentId;
          return (
            <li key={s.id}>
              <button
                ref={cur ? currentRef : undefined}
                type="button"
                aria-current={cur ? "true" : undefined}
                onClick={() => onPick(s.id)}
                className={
                  "flex w-full items-center justify-between gap-2 px-3 py-1 text-left text-body " +
                  (cur ? "bg-row-cursor font-semibold" : "hover:bg-paper-sunk")
                }
              >
                <span className="truncate">{s.name}</span>
                <StatusMark s={s} />
              </button>
            </li>
          );
        })}
        {samples.length === 0 && <li className="px-3 py-2 text-caption">No samples match.</li>}
      </ul>
    </>
  );
}
