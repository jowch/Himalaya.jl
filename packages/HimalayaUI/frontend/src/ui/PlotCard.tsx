import type { ReactNode } from "react";

export function PlotCard({ children, notice }: { children: ReactNode; notice?: ReactNode }): JSX.Element {
  return (
    <div className="relative overflow-hidden rounded-sm border border-hair bg-plate">
      {children}
      {notice && (
        <p role="status" data-testid="workspace-plot-notice" className="pointer-events-none absolute inset-x-0 top-3 mx-auto w-fit max-w-[80%] rounded-sm border border-hair bg-paper px-2.5 py-1 text-center text-caption text-ink-soft">
          {notice}
        </p>
      )}
    </div>
  );
}
