import type { ReactNode } from "react";

export function PlotCard({ children, notice, corner, legend }: { children: ReactNode; notice?: ReactNode; corner?: ReactNode; legend?: ReactNode }): JSX.Element {
  return (
    <div className="relative rounded-sm border border-hair bg-plate">
      {legend && <div data-testid="plot-card-legend" className="relative z-20 px-3 pt-2">{legend}</div>}
      <div className="relative overflow-hidden rounded-b-sm">
        {children}
        {corner && <div data-testid="plot-card-corner" className="absolute right-3 top-6 z-10 flex opacity-40 transition-opacity duration-150 hover:opacity-100 has-focus-visible:opacity-100">{corner}</div>}
        {notice && (
          <p role="status" data-testid="workspace-plot-notice" className={"pointer-events-none absolute inset-x-0 mx-auto w-fit max-w-[80%] rounded-sm border border-hair bg-paper px-2.5 py-1 text-center text-caption text-ink-soft " + (corner ? "top-12" : "top-3")}>
            {notice}
          </p>
        )}
      </div>
    </div>
  );
}
