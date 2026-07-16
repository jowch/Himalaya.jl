import type { ReactNode } from "react";
import { RailSection } from "./RailSection";
import { cx } from "../../lib/cx";


export interface AssignmentRailProps {
  /** The AssignmentCart element (with its PhaseBlock children + onCustomIndex). */
  assignment: ReactNode;
  assignmentCount?: ReactNode;
  /** Optional status note under the cart (mockup .rail-note#assign-note). */
  assignmentNote?: ReactNode;
  /** The CandidateList element (with its CandidateRow children). */
  candidates: ReactNode;
  candidatesNote?: ReactNode;
  className?: string;
}

export function AssignmentRail({
  assignment,
  assignmentCount,
  assignmentNote,
  candidates,
  candidatesNote,
  className,
}: AssignmentRailProps): JSX.Element {
  return (
    <aside
      data-testid="assignment-rail"
      className={cx(
        "bg-paper-sunk border-l border-hair overflow-y-auto p-5 pb-7 flex flex-col gap-5",
        className,
      )}
    >
      <RailSection
        label="Assignment"
        {...(assignmentCount != null ? { count: assignmentCount } : {})}
        {...(assignmentNote != null ? { note: assignmentNote } : {})}
      >
        {assignment}
      </RailSection>
      <RailSection
        label="Candidates"
        {...(candidatesNote != null ? { note: candidatesNote } : {})}
      >
        {candidates}
      </RailSection>
    </aside>
  );
}
