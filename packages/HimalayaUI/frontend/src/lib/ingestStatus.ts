import type { IngestStatus } from "../api";
import type { IngestProgressStatus } from "../state";

export function effectiveIngestStatus(
  inFlight: IngestProgressStatus | undefined,
  persisted: IngestStatus | undefined,
): IngestStatus {
  if (persisted === "complete" || persisted === "failed") return persisted;
  return inFlight ?? persisted ?? "idle";
}