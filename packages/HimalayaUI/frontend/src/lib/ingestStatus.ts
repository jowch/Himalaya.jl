import type { IngestStatus } from "../api";
import type { IngestProgressStatus } from "../state";

/**
 * Reconcile the live SSE overlay (`inFlight`) with the persisted resting truth
 * (`ingest_status`).
 *
 * A TERMINAL persisted state (`complete`/`failed`) WINS over the overlay. The
 * `/api/events` stream is at-most-once: an EventSource reconnect mid-scan (the
 * benign EPIPE disconnects seen in the live walk) drops the `ingest_complete`
 * frame, leaving `ingestInFlight` stuck at `"scanning"` forever — the UI never
 * flips out of the scanning surface until a manual reload re-derives state from
 * the DB. Paired with `useExperiment`'s scoped `refetchInterval` (which polls the
 * resting status while a scan is active), this makes the terminal state self-heal
 * without a reload.
 *
 * Otherwise the overlay leads: it bridges the gap right after Approve, before the
 * experiment row has refetched to `"scanning"`, and carries live progress.
 */
export function effectiveIngestStatus(
  inFlight: IngestProgressStatus | undefined,
  persisted: IngestStatus | undefined,
): IngestStatus {
  // A RESCAN ("analyzing") is driven entirely by the overlay: the persisted row
  // stays "complete" throughout, so the overlay is the source of truth and must
  // win even over a terminal persisted state.
  if (inFlight === "analyzing") return "analyzing";
  // An INITIAL scan transitions the persisted row scanning→complete/failed, so a
  // TERMINAL persisted state is authoritative — it overrides a stale "scanning"
  // overlay left behind by a dropped `ingest_complete` frame (8c).
  if (persisted === "complete" || persisted === "failed") return persisted;
  return inFlight ?? persisted ?? "idle";
}
