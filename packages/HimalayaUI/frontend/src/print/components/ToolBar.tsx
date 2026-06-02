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
      className={`flex items-center gap-2 flex-shrink-0${className ? ` ${className}` : ""}`}
    >
      {children}
    </div>
  );
}
