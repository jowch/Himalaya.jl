import type { ReactNode } from "react";

export interface ToolBarProps {
  children: ReactNode;
  /** PLACEMENT-ONLY (margin / flex position within the plate head). */
  className?: string;
}

export function ToolBar({ children, className }: ToolBarProps): JSX.Element {
  return (
    <div
      role="toolbar"
      // flex-wrap + justify-end: when the head gets tight (~1120px) WHOLE chips
      // drop to a second row, right-aligned, instead of each button's LABEL
      // wrapping to two lines (which left an uneven control row). whitespace-nowrap
      // on the Button / segment primitives keeps each label single-line.
      className={`flex flex-wrap items-center justify-end gap-2 flex-shrink-0${className ? ` ${className}` : ""}`}
    >
      {children}
    </div>
  );
}
