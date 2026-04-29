import { useEffect, useMemo, useRef, useState } from "react";
import { useAppState } from "../state";
import { useIndices, usePeaks, useExposures, useSamples } from "../queries";
import { phaseColor } from "../phases";
import type { IndexEntry, Peak, Exposure, Sample } from "../api";

type PickerRow =
  | { kind: "index";    item: IndexEntry }
  | { kind: "peak";     item: Peak }
  | { kind: "exposure"; item: Exposure }
  | { kind: "sample";   item: Sample };

interface MentionPickerProps {
  query: string;
  onSelect: (token: string) => void;
  onDismiss: () => void;
}

function rowToken(row: PickerRow): string {
  switch (row.kind) {
    case "index":    return `[[index:${row.item.id}]]`;
    case "peak":     return `[[peak:${row.item.id}]]`;
    case "exposure": return `[[exposure:${row.item.id}]]`;
    case "sample":   return `[[sample:${row.item.id}]]`;
  }
}

function matchesQuery(text: string, q: string): boolean {
  return text.toLowerCase().includes(q.toLowerCase());
}

function parseQuery(raw: string): { type: string | null; rest: string } {
  const m = raw.match(/^(\w+):(.*)$/);
  if (m) return { type: m[1]!.toLowerCase(), rest: m[2] ?? "" };
  return { type: null, rest: raw };
}

function rowLabel(row: PickerRow): string {
  switch (row.kind) {
    case "index":    return `${row.item.phase} · ${(row.item.score ?? 0).toFixed(2)}`;
    case "peak":     return `q = ${row.item.q.toFixed(3)}`;
    case "exposure": return row.item.filename ?? `exposure ${row.item.id}`;
    case "sample":   return row.item.name ?? row.item.label ?? `sample ${row.item.id}`;
  }
}

function rowMeta(row: PickerRow): string | null {
  switch (row.kind) {
    case "index":    return `score ${(row.item.score ?? 0).toFixed(2)}`;
    case "peak":     return `${row.item.source} · prom ${(row.item.prominence ?? 0).toFixed(1)}`;
    case "exposure": return row.item.status ?? null;
    case "sample":   return null;
  }
}

export function MentionPicker({ query, onSelect, onDismiss }: MentionPickerProps): JSX.Element {
  const activeExposureId   = useAppState((s) => s.activeExposureId);
  const activeSampleId     = useAppState((s) => s.activeSampleId);
  const activeExperimentId = useAppState((s) => s.activeExperimentId);

  const indicesQ   = useIndices(activeExposureId);
  const peaksQ     = usePeaks(activeExposureId);
  const exposuresQ = useExposures(activeSampleId);
  const samplesQ   = useSamples(activeExperimentId ?? 0);

  const [activeIdx, setActiveIdx] = useState(0);
  const listRef = useRef<HTMLDivElement>(null);
  const scrollRef = useRef<HTMLDivElement>(null);
  const rowRefs = useRef<(HTMLDivElement | null)[]>([]);

  const { type: typeFilter, rest: searchText } = useMemo(() => parseQuery(query), [query]);

  const rows = useMemo((): PickerRow[] => {
    const q = searchText;
    const all: PickerRow[] = [];

    const wantAll   = typeFilter === null;
    const wantIndex = wantAll || typeFilter === "index";
    const wantPeak  = wantAll || typeFilter === "peak";
    const wantExp   = wantAll || typeFilter === "exposure";
    const wantSamp  = wantAll || typeFilter === "sample";

    if (wantIndex) {
      (indicesQ.data ?? [])
        .filter((ix) => !q ||
          matchesQuery(ix.phase, q) ||
          matchesQuery((ix.score ?? 0).toFixed(2), q) ||
          ix.predicted_q.some((pq) => matchesQuery(pq.toFixed(3), q)))
        .forEach((ix) => all.push({ kind: "index", item: ix }));
    }
    if (wantPeak) {
      (peaksQ.data ?? [])
        .filter((pk) => !q || matchesQuery(pk.q.toFixed(3), q))
        .forEach((pk) => all.push({ kind: "peak", item: pk }));
    }
    if (wantExp) {
      (exposuresQ.data ?? [])
        .filter((ex) => !q || matchesQuery(ex.filename ?? "", q))
        .forEach((ex) => all.push({ kind: "exposure", item: ex }));
    }
    if (wantSamp) {
      (samplesQ.data ?? [])
        .filter((sm) => !q || matchesQuery(sm.name ?? sm.label ?? "", q))
        .forEach((sm) => all.push({ kind: "sample", item: sm }));
    }
    return all;
  }, [typeFilter, searchText, indicesQ.data, peaksQ.data, exposuresQ.data, samplesQ.data]);

  useEffect(() => { setActiveIdx(0); }, [rows]);

  // Keep the active row visible when navigating past the visible window.
  useEffect(() => {
    rowRefs.current[activeIdx]?.scrollIntoView({ block: "nearest" });
  }, [activeIdx]);

  useEffect(() => {
    function onKey(e: KeyboardEvent) {
      if (e.key === "Escape") { onDismiss(); return; }
      if (e.key === "ArrowDown") {
        e.preventDefault();
        setActiveIdx((i) => Math.min(i + 1, rows.length - 1));
      }
      if (e.key === "ArrowUp") {
        e.preventDefault();
        setActiveIdx((i) => Math.max(i - 1, 0));
      }
      if (e.key === "Enter" && rows[activeIdx]) {
        e.preventDefault();
        onSelect(rowToken(rows[activeIdx]!));
      }
    }
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [rows, activeIdx, onSelect, onDismiss]);

  return (
    <div
      ref={listRef}
      role="listbox"
      className="absolute bottom-full left-0 right-0 mb-1 z-20
                 bg-bg-elevated border border-border rounded-lg overflow-hidden shadow-xl"
    >
      <div className="px-3 py-1.5 border-b border-border bg-bg
                      text-xs text-fg-dim flex justify-between">
        <span>@{query || "…"}</span>
        <span className="text-fg-dim opacity-50">↑↓ navigate · Enter select · Esc dismiss</span>
      </div>
      {rows.length === 0 ? (
        <div className="px-3 py-2 text-xs text-fg-dim">No results</div>
      ) : (
        <>
          <div
            ref={scrollRef}
            className="max-h-[176px] overflow-y-auto"
          >
            {rows.map((row, i) => {
              const color = row.kind === "index" ? phaseColor(row.item.phase) : undefined;
              const meta  = rowMeta(row);
              return (
                <div
                  key={`${row.kind}:${row.item.id}`}
                  ref={(el) => { rowRefs.current[i] = el; }}
                  role="option"
                  aria-selected={i === activeIdx}
                  onClick={() => onSelect(rowToken(row))}
                  className={`px-3 py-1.5 cursor-pointer flex justify-between items-center text-sm
                              ${i === activeIdx ? "bg-bg-hover" : "hover:bg-bg-hover"}`}
                >
                  <span style={color ? { color } : undefined}>
                    {rowLabel(row)}
                  </span>
                  {meta && (
                    <span className="text-xs text-fg-dim ml-2">{meta}</span>
                  )}
                </div>
              );
            })}
          </div>
          {rows.length > 5 && (
            <div className="px-3 py-1 border-t border-border bg-bg
                            text-[10px] text-fg-dim opacity-70 text-right">
              {rows.length} matches · scroll or refine query
            </div>
          )}
        </>
      )}
    </div>
  );
}
