import { useEffect } from "react";
import { useParams } from "react-router-dom";
import { useAppState } from "../state";
import { useCorpusSamples } from "../queries";

/**
 * useSyncActiveSampleFromRoute — the I4.1 focus-workspace wiring shim.
 *
 * Maps the `/sample/:sampleId` route parameter onto the Zustand
 * `activeSampleId`, so the reused index components (PlotCard / IndicesCard /
 * PhasePanel — which read `active*Id` directly from the store, not from props)
 * render correctly under the URL-routed focus surface (issue #178 / I4.1).
 *
 * One-way (URL -> store). The reverse direction (store -> URL on in-page
 * sample stepping) is owned by I4.2/I4.4 when the focus surface gains its own
 * navigation affordances; reflecting it here would re-introduce the
 * `useUrlFromState` resolving-handshake the corpus shell was built to avoid.
 *
 * The id is validated against the corpus-samples cache so a bogus deep-link
 * (`/sample/999999`) never seeds a non-existent id. `setActiveSample` is
 * called ONLY when the id actually changes, because that action cascades
 * `activeExposureId: undefined` — calling it every render would clobber the
 * exposure `useAutoPickExposure` just selected.
 */
export function useSyncActiveSampleFromRoute(): void {
  const { sampleId: sampleIdParam } = useParams<{ sampleId: string }>();
  const corpusQ = useCorpusSamples();
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const setActiveSample = useAppState((s) => s.setActiveSample);

  useEffect(() => {
    const parsed = Number(sampleIdParam);
    if (!Number.isFinite(parsed)) return; // non-numeric param: ignore
    const samples = corpusQ.data;
    if (samples === undefined) return; // cache not ready yet; re-runs on resolve
    if (!samples.some((s) => s.id === parsed)) return; // unknown sample: ignore
    if (parsed === activeSampleId) return; // no-op: avoid the exposure-clobber cascade
    setActiveSample(parsed);
  }, [sampleIdParam, corpusQ.data, activeSampleId, setActiveSample]);
}
