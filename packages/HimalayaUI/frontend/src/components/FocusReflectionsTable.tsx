import { Skeleton } from "boneyard-js/react";
import { useAppState } from "../state";
import { usePeaks, useIndices, useGroups } from "../queries";
import { phaseColor } from "../phases";
import { HintText } from "./ui";
import type { IndexEntry, Peak } from "../api";

/**
 * FocusReflectionsTable — the third surface of the focus-workspace q-link
 * triple (#209, completes #180). One row per peak; columns are
 *
 *   phase · order · q (Å⁻¹) · d (Å)
 *
 * Row hover sets `hoveredQ` (source), and a row whose q matches the global
 * `hoveredQ` is the sink (data-hot="true"). This wires through the SAME
 * channel that #180 shipped between the trace peaks and the detector rings —
 * no new state field, no new action. Reads `hoveredQ` from the Zustand store
 * and pushes through the same `setHoveredQ` action the trace and the rings
 * already use.
 *
 * Sits in the lower row beside FocusDetectorPanel — `focus-workspace.html`'s
 * `.lower` grid puts the detector square on the left and the reflections
 * table on the right, both keyed to the detector height.
 *
 * Naming note: SAXS peaks do not carry an (h,k,l) Miller triple in our data
 * model — peaks are 1-based positions within the phase's defining ratio
 * series (see `lib/seriesRatio.ts` and `src/phase.jl phaseratios`). The
 * mockup column header is "hkl" but the underlying data is a position. We
 * surface the position under the column header "order" so the rendered
 * value matches the user's understanding (1st reflection, 2nd reflection…).
 */
export function FocusReflectionsTable(): JSX.Element {
  const activeExposureId = useAppState((s) => s.activeExposureId);
  const hoveredQ         = useAppState((s) => s.hoveredQ);
  const setHoveredQ      = useAppState((s) => s.setHoveredQ);

  const peaksQ   = usePeaks(activeExposureId);
  const indicesQ = useIndices(activeExposureId);
  const groupsQ  = useGroups(activeExposureId);

  const peaks   = peaksQ.data ?? [];
  const indices = indicesQ.data ?? [];
  const groups  = groupsQ.data ?? [];

  const activeGroup = groups.find((g) => g.active);
  const memberIds = new Set(activeGroup?.members ?? []);
  const activeIndices = indices.filter((ix) => memberIds.has(ix.id));

  // For each peak, the active index that CLAIMS it (independent phases
  // coexist; same-peak claims are mutually exclusive per `auto_group` /
  // `remove_subsets`, so the first hit is the only hit).
  const claimOf = (peak: Peak): { index: IndexEntry; position: number } | null => {
    for (const ix of activeIndices) {
      const hit = ix.peaks.find((p) => p.peak_id === peak.id);
      if (hit) return { index: ix, position: hit.ratio_position };
    }
    return null;
  };

  // Tolerance for matching hoveredQ to a peak's q. The q-link channel has two
  // tolerance formulas in flight today, not one shared constant:
  //   - the ring overlay uses span-relative — 2% of (qHi - qLo), clamped at a
  //     small floor (DetectorRingOverlay.tsx:73-74 + lib/qRing nearestRingQ);
  //   - the trace peaks use per-peak relative — peak.q * Q_LINK_REL_TOL where
  //     Q_LINK_REL_TOL = 0.01 (TraceViewer.tsx:51, 450-451).
  // The table mirrors the ring formula on purpose: this matches the table's
  // reading order (low-q first, like the rings ordered by radius) and gives a
  // uniform hover band across the row list, where the trace's per-peak rule
  // would shrink the hit-region near q ≈ 0. Both formulas absorb float noise
  // on a hoveredQ that arrives carrying an exact peak.q from a sibling
  // surface; the two are not unified, and the call site that lands on which
  // formula is a current quirk of the codebase, not a stated invariant.
  const peakQs = peaks.map((p) => p.q);
  const qLo = peakQs.length ? Math.min(...peakQs) : 0;
  const qHi = peakQs.length ? Math.max(...peakQs) : 1;
  const tol = Math.max((qHi - qLo) * 0.02, 1e-6);

  // Peaks rendered low-q first, the mockup's reading order. (`peaks` order
  // from the backend is by id; sort by q for the table.)
  const rows = [...peaks].sort((a, b) => a.q - b.q);

  // Footer summary — mirrors mockup `.refl-foot`: "N of M peaks indexed
  // across K phases".
  const covered = peaks.filter((p) => claimOf(p) != null).length;
  const phaseCount = activeIndices.length;

  const body = (() => {
    if (activeExposureId === undefined) {
      return (
        <div className="flex-1 flex items-center justify-center">
          <HintText>Pick a sample to see its reflections.</HintText>
        </div>
      );
    }
    if (peaks.length === 0) {
      return (
        <div className="flex-1 flex items-center justify-center">
          <HintText>No peaks for this exposure.</HintText>
        </div>
      );
    }
    return (
      <>
        <div className="flex-1 min-h-0 overflow-y-auto">
          <table className="w-full border-collapse" data-testid="focus-reflections">
            <thead>
              <tr>
                <th className="sticky top-0 bg-plate border-b border-hair-strong
                               px-2 pb-1.5 pt-1 text-left text-[9.5px] font-bold
                               uppercase tracking-[0.07em] text-ink-faint">
                  phase
                </th>
                <th className="sticky top-0 bg-plate border-b border-hair-strong
                               px-2 pb-1.5 pt-1 text-left text-[9.5px] font-bold
                               uppercase tracking-[0.07em] text-ink-faint">
                  order
                </th>
                <th className="sticky top-0 bg-plate border-b border-hair-strong
                               px-2 pb-1.5 pt-1 text-right text-[9.5px] font-bold
                               uppercase tracking-[0.07em] text-ink-faint font-mono">
                  q (Å⁻¹)
                </th>
                <th className="sticky top-0 bg-plate border-b border-hair-strong
                               px-2 pb-1.5 pt-1 text-right text-[9.5px] font-bold
                               uppercase tracking-[0.07em] text-ink-faint font-mono">
                  d (Å)
                </th>
              </tr>
            </thead>
            <tbody>
              {rows.map((peak) => {
                const claim = claimOf(peak);
                const indexed = claim != null;
                const color = indexed ? phaseColor(claim.index.phase) : "var(--color-ink-faint)";
                // q-link sink: row whose q matches the ephemeral hoveredQ.
                // (Source uses a guarded leave like DetectorRingOverlay so a
                // ring → row hand-off doesn't clobber.)
                const hot = hoveredQ !== undefined
                  && Math.abs(peak.q - hoveredQ) <= tol;
                const dSpacing = peak.q > 0 ? (2 * Math.PI) / peak.q : null;
                return (
                  <tr
                    key={peak.id}
                    data-testid={`reflection-row-${peak.id}`}
                    data-peak-q={peak.q}
                    data-hot={hot ? "true" : "false"}
                    className="cursor-pointer border-b border-hair last:border-b-0
                               transition-colors"
                    style={{
                      background: hot
                        ? "color-mix(in oklab, var(--color-accent) 9%, transparent)"
                        : "transparent",
                    }}
                    onMouseEnter={() => setHoveredQ(peak.q)}
                    onMouseLeave={() =>
                      setHoveredQ(
                        useAppState.getState().hoveredQ === peak.q
                          ? undefined
                          : useAppState.getState().hoveredQ,
                      )}
                  >
                    <td className="px-2 py-2 text-xs">
                      <span
                        aria-hidden
                        className="mr-1.5 inline-block h-[9px] w-[9px] rounded-full
                                   align-[-1px]"
                        style={{ background: color }}
                      />
                      <span
                        className="text-[9.5px] font-bold tracking-[0.03em]"
                        style={{ color }}
                      >
                        {indexed ? claim.index.phase : "unindexed"}
                      </span>
                    </td>
                    <td className="px-2 py-2">
                      <span
                        className={
                          indexed
                            ? "font-mono text-xs font-bold text-ink"
                            : "font-mono text-xs text-ink-faint"
                        }
                      >
                        {indexed ? `#${claim.position}` : "—"}
                      </span>
                    </td>
                    <td className="px-2 py-2 text-right font-mono text-xs tabular-nums
                                   text-ink">
                      {peak.q.toFixed(4)}
                    </td>
                    <td className="px-2 py-2 text-right font-mono text-xs tabular-nums
                                   text-ink">
                      {dSpacing != null ? dSpacing.toFixed(1) : "—"}
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
        <div
          data-testid="focus-reflections-foot"
          className="mt-[11px] pt-[10px] border-t border-hair text-[11px] text-ink-faint
                     shrink-0"
        >
          <span className="font-semibold text-ink-soft">{covered} of {peaks.length}</span>
          {" peaks indexed across "}
          <span className="font-semibold text-ink-soft">
            {phaseCount} phase{phaseCount === 1 ? "" : "s"}
          </span>
          .
        </div>
      </>
    );
  })();

  return (
    // R3-N3 sibling: the reflections panel is plate-like and lifts with the
    // detector panel on its left. Same `.card` Plate Lift treatment.
    <section
      data-testid="focus-reflections-panel"
      className="card flex min-h-0 flex-col p-4"
    >
      <div className="card-header flex items-center justify-between gap-3">
        <span className="text-meta uppercase tracking-wider">Reflections</span>
      </div>
      <Skeleton
        name="focus-reflections"
        className="flex-1 min-h-0 flex flex-col"
        loading={
          activeExposureId !== undefined
          && (peaksQ.isLoading || indicesQ.isLoading || groupsQ.isLoading)
        }
        stagger={50}
        transition={200}
        fallback={
          <div className="flex-1 flex items-center justify-center">
            <HintText>Loading reflections…</HintText>
          </div>
        }
      >
        {body}
      </Skeleton>
    </section>
  );
}
