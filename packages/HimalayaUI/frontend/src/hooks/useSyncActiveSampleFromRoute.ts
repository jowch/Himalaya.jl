import { useEffect } from "react";
import { useParams } from "react-router-dom";
import { useAppState } from "../state";
import { useCorpusSamples } from "../queries";
import type { CorpusSample } from "../api";

/**
 * Route-resolution status for the `/sample/:sampleId` param.
 * "pending": the corpus cache has not resolved yet, so the param cannot be
 * judged. "found": the param names a known sample. "unknown": the cache is
 * ready and the param is non-numeric, or numeric but matching no sample.
 */
export type RouteSampleStatus = "pending" | "found" | "unknown";

/**
 * The single shared predicate behind RouteSampleStatus: the exact parse plus
 * cache lookup the seeding hook uses. Pure, so honest surfaces rendered
 * outside the routed element (the sample stepper in the app shell) can judge the
 * same param against the same cache without duplicating the logic.
 */
export function resolveRouteSampleStatus(
  sampleIdParam: string | undefined,
  samples: CorpusSample[] | undefined,
): RouteSampleStatus {
  if (samples === undefined) return "pending";
  const parsed = Number(sampleIdParam);
  if (!Number.isFinite(parsed)) return "unknown";
  return samples.some((s) => s.id === parsed) ? "found" : "unknown";
}

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
 *
 * Returns the RouteSampleStatus from `resolveRouteSampleStatus`, the same
 * predicate the seeding effect mirrors. A bogus param never seeds the store,
 * but mid-session the previous `activeSampleId` survives; the caller
 * (FocusPage) uses the "unknown" status to render not-found instead of
 * silently showing the previous sample under the wrong URL.
 */
export function useSyncActiveSampleFromRoute(): RouteSampleStatus {
  const { sampleId: sampleIdParam } = useParams<{ sampleId: string }>();
  const corpusQ = useCorpusSamples();
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const setActiveSample = useAppState((s) => s.setActiveSample);

  // Single parse + cache handle: the seeding effect and the returned status
  // both derive from these, so the two can never diverge.
  const parsed = Number(sampleIdParam);
  const samples = corpusQ.data;

  useEffect(() => {
    if (!Number.isFinite(parsed)) return; // non-numeric param: never seed
    if (samples === undefined) return; // cache not ready yet; re-runs on resolve
    if (!samples.some((s) => s.id === parsed)) return; // unknown sample: never seed
    if (parsed === activeSampleId) return; // no-op: avoid the exposure-clobber cascade
    setActiveSample(parsed);
  }, [parsed, samples, activeSampleId, setActiveSample]);

  return resolveRouteSampleStatus(sampleIdParam, samples);
}
