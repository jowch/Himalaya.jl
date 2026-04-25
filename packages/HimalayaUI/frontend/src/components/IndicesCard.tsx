import { useMemo } from "react";
import { PhasePanel } from "./PhasePanel";
import { MillerPlot } from "./MillerPlot";
import { useAppState } from "../state";
import { useIndices, useGroups } from "../queries";

/**
 * IndicesCard — right card on the Index page. Top: PhasePanel (Active set +
 * Candidates). Bottom: a small MillerPlot showing √N vs q for the active
 * indices (a "did the indexing land on a clean sqrt-of-integers line?" sanity
 * check). The Miller plot used to live as an inset over the trace itself, but
 * any in-plot position covers something when the user zooms — moving it
 * outside the trace plot eliminates that conflict entirely.
 */
export function IndicesCard(): JSX.Element {
  const exposureId    = useAppState((s) => s.activeExposureId);
  const hoveredIndexId = useAppState((s) => s.hoveredIndexId);
  const indicesQ      = useIndices(exposureId);
  const groupsQ       = useGroups(exposureId);

  const { activeGroupIndices, hoveredIndex } = useMemo(() => {
    const indices = indicesQ.data ?? [];
    const active  = (groupsQ.data ?? []).find((g) => g.active);
    const activeGroupIndices = (active?.members ?? [])
      .map((id) => indices.find((i) => i.id === id))
      .filter((i): i is NonNullable<typeof i> => i != null);
    const hoveredIndex = hoveredIndexId != null
      ? indices.find((i) => i.id === hoveredIndexId)
      : undefined;
    return { activeGroupIndices, hoveredIndex };
  }, [indicesQ.data, groupsQ.data, hoveredIndexId]);

  return (
    <div data-testid="indices-card" className="h-full min-h-0 flex flex-col">
      <div className="flex-1 min-h-0">
        <PhasePanel exposureId={exposureId} />
      </div>
      <div
        data-testid="miller-panel"
        className="shrink-0 border-t border-border bg-bg-elevated/60"
      >
        <div className="h-[18px] px-3 flex items-center justify-between
                        text-[9.5px] uppercase tracking-wider text-fg-dim
                        border-b border-border">
          <span>√N · q sanity</span>
        </div>
        <div className="h-[140px]">
          <MillerPlot indices={activeGroupIndices} hoveredIndex={hoveredIndex} />
        </div>
      </div>
    </div>
  );
}
