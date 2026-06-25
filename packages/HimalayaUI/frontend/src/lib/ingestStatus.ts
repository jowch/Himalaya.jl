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
  // A TERMINAL persisted state (complete/failed) is authoritative and overrides a
  // stale overlay. BOTH an initial scan and a rescan transition the persisted row
  // scanning→complete/failed (routes_experiments.jl), so a dropped terminal SSE
  // frame can strand the overlay at "scanning" (initial) OR "analyzing" (rescan);
  // either way the polled terminal row wins and self-heals the surface.
  if (persisted === "complete" || persisted === "failed") return persisted;
  // Live scan/rescan: the overlay carries the phase distinction (scanning vs
  // analyzing) and leads; persisted "scanning" bridges the pre-overlay gap.
  return inFlight ?? persisted ?? "idle";
}
