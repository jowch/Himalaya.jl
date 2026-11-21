import { useEffect, useRef, type ReactNode } from "react";
import type { CorpusSample } from "../../api";
import { sampleStatus } from "../model";

interface Props {
  samples: readonly CorpusSample[];
  currentId: number | undefined;
  filter: string;
  onFilter: (v: string) => void;
  onPick: (id: number) => void;
  bridgeCounts?: ReadonlyMap<number, number>;
  notice?: ReactNode;
}

function StatusMark({ s }: { s: CorpusSample }): JSX.Element {
  const st = sampleStatus(s);
  if (st.kind === "phase") return <span data-testid="workspace-status" className="text-xs font-semibold text-ink">{st.phase}</span>;
  if (st.kind === "form_factor") return <span data-testid="workspace-status" className="text-xs text-ink-faint">Form factor</span>;
  return <span data-testid="workspace-status" role="img" aria-label="Unreviewed" title="Unreviewed" className="text-xs text-ink-faint">○</span>;
}

export function SampleList({ samples, currentId, filter, onFilter, onPick, bridgeCounts, notice }: Props): JSX.Element {
  const currentRef = useRef<HTMLButtonElement | null>(null);
  useEffect(() => {
    currentRef.current?.scrollIntoView({ block: "nearest" });
  }, [currentId]);
  return (
    <>
      <div className="border-b border-hair px-2.5 pb-2 pt-2.5">
        <input
          id="workspace-filter"
          type="search"
          value={filter}
          onChange={(e) => onFilter(e.target.value)}
          onKeyDown={(e) => {
            if (e.nativeEvent.isComposing) return;
            if (e.key === "Escape") {
              e.preventDefault();
              e.currentTarget.blur();
            } else if (e.key === "Enter") {
              e.preventDefault();
              const first = samples[0];
              if (first) {
                e.currentTarget.blur();
                onPick(first.id);
              }
            }
          }}
          placeholder="Filter  /"
          aria-label="Filter samples"
          className="w-full rounded-sm border border-hair bg-plate px-2 py-[5px] text-base text-ink placeholder:text-ink-faint focus:outline-2 focus:-outline-offset-1 focus:outline-accent"
        />
      </div>
      {notice}
      <ul className="min-h-0 flex-1 overflow-y-auto" data-testid="workspace-samples">
        {samples.map((s) => {
          const cur = s.id === currentId;
          const bridges = bridgeCounts?.get(s.id) ?? 0;
          return (
            <li key={s.id}>
              <button
                ref={cur ? currentRef : undefined}
                type="button"
                aria-current={cur ? "true" : undefined}
                data-sample-id={s.id}
                onClick={() => onPick(s.id)}
                className={
                  "flex w-full items-center justify-between gap-2 px-2.5 py-[5px] text-left text-base font-semibold text-ink " +
                  (cur ? "bg-row-cursor shadow-[inset_3px_0_0_var(--color-accent)]" : "hover:bg-paper-sunk")
                }
              >
                <span className="flex min-w-0 items-center gap-2">
                  <span className="truncate">{s.name}</span>
                  {bridges > 0 && (
                    <span data-bridged title={bridges === 1 ? "1 bridge" : `${bridges} bridges`} className="shrink-0 text-xs text-ink-faint tabular-nums">
                      <span aria-hidden="true">⇄{bridges}</span>
                      <span className="sr-only">{bridges === 1 ? "1 bridge" : `${bridges} bridges`}</span>
                    </span>
                  )}
                </span>
                <StatusMark s={s} />
              </button>
            </li>
          );
        })}
        {samples.length === 0 && filter.trim() !== "" && <li className="px-3 py-2 text-caption">No samples match.</li>}
      </ul>
    </>
  );
}
