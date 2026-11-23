import { useEffect, useRef, useState, type FocusEvent, type ReactNode } from "react";

export const NAV_PEEK_DELAY_MS = 250;
export const NAV_PEEK_LEAVE_MS = 200;

interface Props {
  header: ReactNode;
  nav: ReactNode;
  main: ReactNode;
  tools: ReactNode;
  theme?: string;
  overlay?: ReactNode;
  navHidden?: boolean;
  onShowNav?: () => void;
}

export function WorkspaceLayout({ header, nav, main, tools, theme = "system", overlay, navHidden = false, onShowNav }: Props): JSX.Element {
  const [hovered, setHovered] = useState(false);
  const [focused, setFocused] = useState(false);
  const timer = useRef<number | undefined>(undefined);
  const clear = (): void => {
    if (timer.current !== undefined) window.clearTimeout(timer.current);
    timer.current = undefined;
  };
  const later = (open: boolean, ms: number): void => {
    clear();
    timer.current = window.setTimeout(() => {
      timer.current = undefined;
      setHovered(open);
    }, ms);
  };
  useEffect(() => clear, []);
  useEffect(() => {
    clear();
    setHovered(false);
  }, [navHidden]);
  const peek = navHidden && (hovered || focused);
  const onBlur = (e: FocusEvent<HTMLElement>): void => {
    if (!(e.relatedTarget instanceof Node && e.currentTarget.contains(e.relatedTarget))) setFocused(false);
  };
  return (
    <div className="grid h-full grid-rows-[auto_1fr] bg-paper text-ink" data-testid="workspace" data-theme={theme}>
      <header className="flex min-w-0 items-center gap-3.5 border-b border-hair bg-plate px-3.5 py-2">{header}</header>
      <div className={"relative grid min-h-0 overflow-clip " + (navHidden ? "grid-cols-[14px_minmax(0,1fr)_300px]" : "grid-cols-[270px_minmax(0,1fr)_300px]")}>
        {navHidden && (
          <button
            type="button"
            aria-label="Show the sample list"
            data-testid="sample-list-handle"
            onClick={onShowNav}
            onMouseEnter={() => later(true, NAV_PEEK_DELAY_MS)}
            onMouseLeave={() => later(false, NAV_PEEK_LEAVE_MS)}
            className="flex min-h-0 justify-center border-r border-hair bg-paper pt-2.5 text-sm text-ink-faint hover:bg-paper-sunk hover:text-ink focus-visible:outline-2 focus-visible:-outline-offset-2 focus-visible:outline-accent"
          >
            <span aria-hidden="true">›</span>
          </button>
        )}
        <nav
          aria-label="Samples"
          data-peek={peek ? "true" : undefined}
          onMouseEnter={navHidden ? clear : undefined}
          onMouseLeave={navHidden ? () => later(false, NAV_PEEK_LEAVE_MS) : undefined}
          onFocus={navHidden ? () => setFocused(true) : undefined}
          onBlur={navHidden ? onBlur : undefined}
          className={
            "flex min-h-0 flex-col overflow-hidden border-r border-hair " +
            (navHidden
              ? "absolute inset-y-0 left-0 z-30 w-[270px] bg-paper transition-transform duration-150 " + (peek ? "translate-x-0 shadow-lg" : "-translate-x-full")
              : "")
          }
        >
          {nav}
        </nav>
        <main className="flex min-h-0 min-w-0 flex-col overflow-y-auto">{main}</main>
        <aside aria-label="Tools" className="flex min-h-0 flex-col overflow-y-auto border-l border-hair">{tools}</aside>
      </div>
      {overlay}
    </div>
  );
}
