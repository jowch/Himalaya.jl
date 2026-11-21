import type { ReactNode } from "react";

export function PlotCard({ children }: { children: ReactNode }): JSX.Element {
  return <div className="overflow-hidden rounded-sm border border-hair bg-plate">{children}</div>;
}
