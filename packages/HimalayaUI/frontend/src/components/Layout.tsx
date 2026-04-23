import type { ReactNode } from "react";

export interface LayoutProps {
  left: ReactNode;
  centerTop?: ReactNode;
  centerBottom?: ReactNode;
  rightTop?: ReactNode;
  rightBottom?: ReactNode;
}

function Placeholder({ label }: { label: string }): JSX.Element {
  return <p className="muted placeholder">{label}</p>;
}

export function Layout({
  left, centerTop, centerBottom, rightTop, rightBottom,
}: LayoutProps): JSX.Element {
  return (
    <main className="layout">
      <aside className="col col-left">{left}</aside>
      <div className="col col-center">
        <section className="pane pane-center-top">
          {centerTop ?? <Placeholder label="Trace viewer (Plan 4)" />}
        </section>
        <section className="pane pane-center-bottom">
          {centerBottom ?? <Placeholder label="Properties panel (Plan 6)" />}
        </section>
      </div>
      <div className="col col-right">
        <section className="pane pane-right-top">
          {rightTop ?? <Placeholder label="Miller-index plot (Plan 5)" />}
        </section>
        <section className="pane pane-right-bottom">
          {rightBottom ?? <Placeholder label="Phase panel (Plan 5)" />}
        </section>
      </div>
    </main>
  );
}
