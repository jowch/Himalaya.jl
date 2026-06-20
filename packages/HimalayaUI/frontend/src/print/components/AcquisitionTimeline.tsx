import type { JSX } from "react";
import { AcquisitionChart } from "../plot/AcquisitionChart"; // render layer (appearance-exempt)

export interface AcqSession { label: string; loadFrameCounts: number[] }

export interface AcquisitionTimelineProps {
  sessions: AcqSession[];
  className?: string;
}

/** Per-session cluster of per-load exposure-count bars, drawn on the shared
 *  trace-plot family. Placement wrapper; appearance lives in the exempt render
 *  layer (print/plot/AcquisitionChart). */
export function AcquisitionTimeline({ sessions, className }: AcquisitionTimelineProps): JSX.Element {
  return (
    <div className={className} data-testid="acquisition-timeline">
      <AcquisitionChart sessions={sessions} />
    </div>
  );
}
