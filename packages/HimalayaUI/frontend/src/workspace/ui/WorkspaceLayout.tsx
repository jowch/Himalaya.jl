import type { ReactNode } from "react";

interface Props {
  header: ReactNode;
  nav: ReactNode;
  main: ReactNode;
  tools: ReactNode;
  theme?: string;
  overlay?: ReactNode;
}

export function WorkspaceLayout({ header, nav, main, tools, theme = "system", overlay }: Props): JSX.Element {
  return (
    <div className="grid h-full grid-rows-[auto_1fr] bg-paper text-ink" data-testid="workspace" data-workspace-theme={theme}>
      <header className="flex min-w-0 items-center gap-3.5 border-b border-hair bg-plate px-3.5 py-2">{header}</header>
      <div className="grid min-h-0 grid-cols-[270px_minmax(0,1fr)_300px]">
        <nav aria-label="Samples" className="flex min-h-0 flex-col overflow-hidden border-r border-hair">{nav}</nav>
        <main className="flex min-h-0 min-w-0 flex-col overflow-y-auto">{main}</main>
        <aside aria-label="Tools" className="flex min-h-0 flex-col overflow-y-auto border-l border-hair">{tools}</aside>
      </div>
      {overlay}
    </div>
  );
}
