import { useExperiment } from "../queries";
import { useAppState } from "../state";
import { effectiveIngestStatus } from "../lib/ingestStatus";

export interface ScanState {
  loaded: boolean;
  scanning: boolean;
  rescanning: boolean;
  failed: boolean;
}

export function useScanState(experimentId: number): ScanState {
  const inFlight = useAppState((s) => s.ingestInFlight?.[experimentId]);
  const exp = useExperiment(experimentId);
  const eff = effectiveIngestStatus(inFlight?.status, exp.data?.ingest_status);
  const processing = eff === "scanning" || eff === "analyzing";
  const initial = eff === "scanning" && exp.data?.last_scanned_at == null;
  return {
    loaded: exp.data !== undefined,
    scanning: processing && initial,
    rescanning: processing && !initial,
    failed: !processing && eff === "failed",
  };
}
