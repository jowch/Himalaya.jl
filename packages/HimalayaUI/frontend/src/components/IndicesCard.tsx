import { PhasePanel } from "./PhasePanel";
import { useAppState } from "../state";

/**
 * IndicesCard — the focus-workspace rail's "output" surface. It renders the
 * PhasePanel (the "Phase call" block + multi-select "Candidate indexings").
 *
 * The legacy √N·q MillerPlot inset was dropped from the Print rail (R4, finding
 * L-12): it was an extra dark-token surface the mockup rail does not have. The
 * √N-series sanity check now lives in the Phase-call block's series ratio
 * (√2 : √3 : √4) and on each candidate row. The MillerPlot component itself is
 * retained for other consumers; it just no longer hangs under the focus rail.
 */
export function IndicesCard(): JSX.Element {
  const exposureId = useAppState((s) => s.activeExposureId);

  return (
    <div data-testid="indices-card" className="h-full min-h-0 flex flex-col">
      <PhasePanel exposureId={exposureId} />
    </div>
  );
}
