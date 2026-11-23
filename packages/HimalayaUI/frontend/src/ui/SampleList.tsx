import { Fragment, useEffect, useMemo, useRef, useState, type ReactNode } from "react";
import type { CorpusSample } from "../api";
import { sampleStatus } from "../workspace/model";
import { COMPARE_DASH } from "./WorkspacePlot";
import { SELECTED_ROW } from "./selectedRow";
import { Tooltip } from "./Tooltip";
import { isShortened, PhaseShortName } from "./PhaseName";
import { useListHighlight } from "./useListHighlight";
import { Dot } from "./Dot";
import { IconButton } from "./IconButton";
import { releaseMouseFocus, useRovingFocus } from "../hooks/useRovingFocus";

interface Props {
  samples: readonly CorpusSample[];
  setupHeaders?: ReadonlyMap<number, string> | undefined;
  currentId: number | undefined;
  filter: string;
  onFilter: (v: string) => void;
  onPick: (id: number) => void;
  onToggleCompare?: (id: number) => void;
  compareColors?: ReadonlyMap<number, string>;
  hiddenCompares?: ReadonlySet<number>;
  bridgeCounts?: ReadonlyMap<number, number>;
  notice?: ReactNode;
  listHidden?: boolean;
  onToggleList?: () => void;
}

const ANALYSIS_GAP = {
  all_dropped: { mark: "All dropped", why: "Every image is dropped, so nothing is analyzed.", label: "Unreviewed. Every image is dropped, so nothing is analyzed." },
  none_selected: { mark: "None selected", why: "No image is selected, so nothing is analyzed.", label: "Unreviewed. No image is selected, so nothing is analyzed." },
} as const;

function StatusMark({ s }: { s: CorpusSample }): JSX.Element {
  const st = sampleStatus(s);
  if (st.kind === "phase") {
    const name = st.phases.join(" + ");
    const mark = (
      <span data-testid="workspace-status" role="img" aria-label={name} className="flex gap-1 text-xs font-semibold text-ink">
        {st.phases.map((p) => <PhaseShortName key={p} phase={p} />)}
      </span>
    );
    return st.phases.some(isShortened) ? <Tooltip content={name}>{mark}</Tooltip> : mark;
  }
  if (st.kind === "form_factor") return <span data-testid="workspace-status" className="text-xs text-ink-faint">Form factor</span>;
  const gap = s.analysis_gap ? ANALYSIS_GAP[s.analysis_gap] : null;
  if (gap) {
    return (
      <Tooltip content={gap.why}>
        <span data-testid="workspace-status" role="img" aria-label={gap.label} className="text-xs text-ink-faint">{gap.mark}</span>
      </Tooltip>
    );
  }
  return (
    <Tooltip content="Unreviewed">
      <span data-testid="workspace-status" role="img" aria-label="Unreviewed" className="text-xs text-ink-faint">○</span>
    </Tooltip>
  );
}

export function SampleList({ samples, setupHeaders, currentId, filter, onFilter, onPick, onToggleCompare, compareColors, hiddenCompares, bridgeCounts, notice, listHidden = false, onToggleList }: Props): JSX.Element {
  const currentRef = useRef<HTMLButtonElement | null>(null);
  useEffect(() => {
    currentRef.current?.scrollIntoView({ block: "nearest" });
  }, [currentId]);
  const filterRef = useRef<HTMLInputElement | null>(null);
  const [filterFocused, setFilterFocused] = useState(false);
  const ids = useMemo(() => samples.map((s) => s.id), [samples]);
  const highlight = useListHighlight({
    count: samples.length,
    keys: ids,
    resetKey: filter,
    onCommit: (i) => {
      const picked = samples[i];
      if (!picked) return;
      filterRef.current?.blur();
      onPick(picked.id);
    },
    onEscape: () => filterRef.current?.blur(),
  });
  const highlightedId = filterFocused ? samples[highlight.at]?.id : undefined;
  useEffect(() => {
    if (highlightedId !== undefined) document.getElementById(`workspace-sample-${highlightedId}`)?.scrollIntoView({ block: "nearest" });
  }, [highlightedId]);
  const listRef = useRovingFocus<HTMLUListElement>(currentId);
  const hasCurrent = samples.some((s) => s.id === currentId);
  return (
    <>
      <div className="flex items-center gap-1.5 border-b border-hair px-2.5 pb-2 pt-2.5">
        <input
          id="workspace-filter"
          type="search"
          value={filter}
          onChange={(e) => onFilter(e.target.value)}
          ref={filterRef}
          onFocus={() => setFilterFocused(true)}
          onBlur={() => setFilterFocused(false)}
          onKeyDown={(e) => {
            highlight.onKeyDown(e);
          }}
          aria-activedescendant={highlightedId !== undefined ? `workspace-sample-${highlightedId}` : undefined}
          placeholder="Filter  /"
          aria-label="Filter samples"
          className="min-w-0 flex-1 rounded-sm border border-hair bg-plate px-2 py-[5px] text-base text-ink placeholder:text-ink-faint focus:outline-2 focus:-outline-offset-1 focus:outline-accent"
        />
        {onToggleList && (
          <IconButton label={listHidden ? "Keep the sample list open" : "Hide the sample list"} data-testid="sample-list-toggle" onClick={onToggleList}>
            <span aria-hidden="true" className="text-base">{listHidden ? "›" : "‹"}</span>
          </IconButton>
        )}
      </div>
      {notice}
      <ul ref={listRef} className="min-h-0 flex-1 overflow-y-auto" data-testid="workspace-samples">
        {samples.map((s, i) => {
          const cur = s.id === currentId;
          const setupId = s.setup_id ?? undefined;
          const header = setupHeaders && setupId !== undefined && (i === 0 || samples[i - 1]!.setup_id !== s.setup_id)
            ? setupHeaders.get(setupId)
            : undefined;
          const bridges = bridgeCounts?.get(s.id) ?? 0;
          const compareColor = compareColors?.get(s.id);
          const beamOff = s.beam_off_images ?? 0;
          const beamOffText = beamOff === 1 ? "1 image likely beam-off" : `${beamOff} images likely beam-off`;
          return (
            <Fragment key={s.id}>
            {header !== undefined && (
              <li data-testid="setup-header" className="sticky top-0 z-[1] border-b border-hair bg-paper px-2.5 pb-1 pt-2 text-xs font-semibold text-ink-soft">
                {header}
              </li>
            )}
            <li>
              <button
                ref={cur ? currentRef : undefined}
                id={`workspace-sample-${s.id}`}
                type="button"
                data-highlighted={s.id === highlightedId ? "" : undefined}
                aria-current={cur ? "true" : undefined}
                tabIndex={cur || (!hasCurrent && i === 0) ? 0 : -1}
                data-roving="row"
                data-sample-id={s.id}
                data-compared={compareColor !== undefined ? "" : undefined}
                onClick={(e) => {
                  releaseMouseFocus(e);
                  if ((e.metaKey || e.ctrlKey) && onToggleCompare) {
                    e.preventDefault();
                    onToggleCompare(s.id);
                    return;
                  }
                  onPick(s.id);
                }}
                className={
                  "flex w-full items-center justify-between gap-2 px-2.5 py-[5px] text-left text-base text-ink " +
                  (cur ? "font-semibold " + SELECTED_ROW : s.id === highlightedId ? "bg-paper-sunk" : "hover:bg-paper-sunk")
                }
              >
                <span className="flex min-w-0 items-center gap-2">
                  <span className="truncate">{s.name}</span>
                  {bridges > 0 && (
                    <Tooltip content={bridges === 1 ? "1 bridge" : `${bridges} bridges`}>
                      <span data-bridged className="shrink-0 text-xs text-ink-faint tabular-nums">
                        <span aria-hidden="true">⇄{bridges}</span>
                        <span className="sr-only">{bridges === 1 ? "1 bridge" : `${bridges} bridges`}</span>
                      </span>
                    </Tooltip>
                  )}
                </span>
                <span className="flex shrink-0 items-center gap-2">
                  {beamOff > 0 && (
                    <Tooltip content={beamOffText}>
                      <span data-warned="" className="flex">
                        <Dot tone="warning" size="xs" label={beamOffText} />
                      </span>
                    </Tooltip>
                  )}
                  {compareColor !== undefined && (
                    <svg width={14} height={4} role="img" aria-label={hiddenCompares?.has(s.id) ? "Compared, hidden" : "Compared"} style={{ color: compareColor }} className={"shrink-0" + (hiddenCompares?.has(s.id) ? " opacity-45" : "")}>
                      <line x1={0} x2={14} y1={2} y2={2} stroke="currentColor" strokeWidth={2} strokeDasharray={COMPARE_DASH} />
                    </svg>
                  )}
                  <StatusMark s={s} />
                </span>
              </button>
            </li>
            </Fragment>
          );
        })}
        {samples.length === 0 && filter.trim() !== "" && <li className="px-3 py-2 text-caption">No samples match.</li>}
      </ul>
    </>
  );
}
