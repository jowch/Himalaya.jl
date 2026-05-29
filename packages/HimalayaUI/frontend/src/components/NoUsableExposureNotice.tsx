import { Link } from "react-router-dom";

export interface NoUsableExposureNoticeProps {
  sampleId: number;
  /** True when the sample has exposures but all are rejected (vs. has none). */
  allRejected: boolean;
}

/**
 * Shown on the focus workspace when the active sample has no usable exposure to
 * index — every exposure is rejected, or the sample has none at all. The M1
 * corpus→focus doors (`ContactSheetRow` / `LoupeSidebar` / `NavModal`) can land
 * a scientist here, so this says *why* the trace is empty and points back to
 * the loupe, where rejection is reversible (the Drop/Restore + Set-representative
 * controls live there). Without it the plot falls to the generic "No trace data
 * available", which reads as a load failure and offers no way forward — a door
 * into a broken room.
 *
 * The single terracotta link is the rationed grease-pencil mark (DESIGN.md): an
 * empty state with one interactive exit, not decoration.
 */
export function NoUsableExposureNotice({
  sampleId,
  allRejected,
}: NoUsableExposureNoticeProps): JSX.Element {
  return (
    <div
      data-testid="no-usable-exposure"
      className="flex-1 flex flex-col items-center justify-center gap-2 px-6 text-center"
    >
      <p className="text-title text-ink">
        {allRejected ? "Every exposure is rejected" : "No exposures yet"}
      </p>
      <p className="max-w-[42ch] text-sm text-ink-faint">
        {allRejected
          ? "This sample has no accepted exposure to index."
          : "This sample has no exposures to index."}
      </p>
      <Link
        to={`/samples/loupe/${sampleId}`}
        data-testid="no-usable-exposure-loupe-link"
        className="text-sm font-medium text-print-accent hover:underline
                   focus-visible:outline focus-visible:outline-2
                   focus-visible:outline-accent focus-visible:outline-offset-2"
      >
        {allRejected ? "Restore an exposure in the loupe" : "Open the loupe"}
      </Link>
    </div>
  );
}
