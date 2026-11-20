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
      <header className="flex items-center gap-3 border-b border-hair bg-plate px-4 py-2">{header}</header>
      <div className="grid min-h-0 grid-cols-[260px_minmax(0,1fr)_300px]">
        <nav aria-label="Samples" className="flex min-h-0 flex-col border-r border-hair">{nav}</nav>
        <main className="flex min-h-0 min-w-0 flex-col gap-4 overflow-y-auto px-5 py-4">{main}</main>
        <aside aria-label="Tools" className="flex min-h-0 flex-col gap-6 overflow-y-auto border-l border-hair px-4 py-4">{tools}</aside>
      </div>
      {overlay}
    </div>
  );
}
