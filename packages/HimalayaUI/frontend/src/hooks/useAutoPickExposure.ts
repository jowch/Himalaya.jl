import { useEffect } from "react";
import { useAppState } from "../state";
import { useExposures } from "../queries";

/**
 * Auto-select the right exposure for the active sample on the Index page:
 *
 *   1. If the user has marked one for indexing on the Inspect page
 *      (`selected` flag set by `PATCH /api/exposures/:id/select`), follow
 *      it. This is the load-bearing case — without it, marking for
 *      indexing has no visible effect on the Index page.
 *   2. Otherwise keep the current `activeExposureId` if it still exists
 *      on this sample.
 *   3. Otherwise fall back to the first non-rejected exposure.
 *
 * Extracted from PlotCard so the regression test for issue #118
 * (`test/sampleSwitchKeypress.test.tsx`) can subscribe to the same effect
 * the production tree runs, eliminating drift between a hand-maintained
 * mirror and the real auto-pick logic.
 */
export function useAutoPickExposure(activeSampleId: number | undefined): void {
  const activeExposureId = useAppState((s) => s.activeExposureId);
  const setActiveExposure = useAppState((s) => s.setActiveExposure);
  const exposuresQ = useExposures(activeSampleId);
  useEffect(() => {
    // Skip rejected exposures. Hook returns the unfiltered list (one cache
    // row shared with Inspect); the Index page is the only consumer that
    // wants acceptable-only.
    const exposures = (exposuresQ.data ?? []).filter((e) => e.status !== "rejected");
    if (exposures.length === 0) return;
    const flagged = exposures.find((e) => e.selected);
    if (flagged) {
      if (activeExposureId !== flagged.id) setActiveExposure(flagged.id);
      return;
    }
    const stillValid = exposures.some((e) => e.id === activeExposureId);
    if (!stillValid) setActiveExposure(exposures[0]!.id);
  }, [exposuresQ.data, activeExposureId, setActiveExposure]);
}
