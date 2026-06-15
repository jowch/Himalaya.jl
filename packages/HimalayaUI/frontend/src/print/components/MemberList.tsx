import type { ReactNode } from "react";
import { SeriesMemberRow } from "./SeriesMemberRow";

export interface MemberDatum {
  key: string;
  phases: string[];
  variableValue: ReactNode;
  dataLine: ReactNode;
}

export interface MemberListProps {
  /** Members in display order (page reverses so top = high variable). */
  members: MemberDatum[];
  /** Controlled hot key, synced with the waterfall hot row. */
  hoveredKey?: string;
  onHoverMember?: (key?: string) => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function MemberList({ members, hoveredKey, onHoverMember, className }: MemberListProps): JSX.Element {
  return (
    <div data-testid="member-list" className={`flex flex-col gap-0.5${className ? ` ${className}` : ""}`}>
      {members.map((m) => (
        <SeriesMemberRow
          key={m.key}
          phases={m.phases}
          variableValue={m.variableValue}
          dataLine={m.dataLine}
          hot={hoveredKey === m.key}
          {...(onHoverMember ? { onHover: () => onHoverMember(m.key), onLeave: () => onHoverMember(undefined) } : {})}
        />
      ))}
    </div>
  );
}
