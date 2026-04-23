import type { ReactNode } from "react";

export interface LayoutProps {
  left: ReactNode;
  centerTop?: ReactNode;
  centerBottom?: ReactNode;
  rightTop?: ReactNode;
  rightBottom?: ReactNode;
}

const pane = "flex-1 overflow-auto p-4 flex flex-col min-h-0";
const paneBorder = "border-t border-border";
const placeholder = "flex items-center justify-center flex-1 text-fg-dim italic";

function Placeholder({ label }: { label: string }): JSX.Element {
  return <p className={placeholder} data-testid="placeholder">{label}</p>;
}

export function Layout({
  left, centerTop, centerBottom, rightTop, rightBottom,
}: LayoutProps): JSX.Element {
  return (
    <main className="grid grid-cols-[280px_1fr_360px] h-[calc(100vh-44px)] overflow-hidden">
      <aside className="flex flex-col overflow-hidden border-r border-border">{left}</aside>
      <div className="flex flex-col overflow-hidden border-r border-border min-w-0">
        <section className={pane}>{centerTop ?? <Placeholder label="Trace viewer (Plan 4)" />}</section>
        <section className={`${pane} ${paneBorder}`}>{centerBottom ?? <Placeholder label="Properties panel (Plan 6)" />}</section>
      </div>
      <div className="flex flex-col overflow-hidden min-w-0">
        <section className={pane}>{rightTop ?? <Placeholder label="Miller-index plot (Plan 5)" />}</section>
        <section className={`${pane} ${paneBorder}`}>{rightBottom ?? <Placeholder label="Phase panel (Plan 5)" />}</section>
      </div>
    </main>
  );
}
