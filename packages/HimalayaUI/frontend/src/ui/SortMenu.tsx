import { useEffect, useRef, useState, type KeyboardEvent } from "react";
import { Tooltip } from "./Tooltip";

export interface SortChoice<T extends string> {
  sort: T;
  name: string;
  hint: string;
}

interface Props<T extends string> {
  what: string;
  current: T;
  choices: readonly SortChoice<T>[];
  onSort: (sort: T) => void;
}

const STEP: Record<string, number> = { ArrowDown: 1, ArrowUp: -1 };

function SortGlyph(): JSX.Element {
  return (
    <svg width="12" height="12" viewBox="0 0 12 12" aria-hidden="true" className="shrink-0">
      <path d="M3 2v8M1 8l2 2 2-2M9 10V2M7 4l2-2 2 2" fill="none" stroke="currentColor" strokeWidth="1.2" strokeLinecap="round" strokeLinejoin="round" />
    </svg>
  );
}

export function SortMenu<T extends string>({ what, current, choices, onSort }: Props<T>): JSX.Element {
  const [open, setOpen] = useState(false);
  const root = useRef<HTMLDivElement | null>(null);
  const button = useRef<HTMLButtonElement | null>(null);
  const shown = choices.find((c) => c.sort === current) ?? choices[0];
  useEffect(() => {
    if (!open) return;
    root.current?.querySelector<HTMLButtonElement>("[aria-checked='true']")?.focus();
    const onDown = (ev: MouseEvent): void => {
      if (root.current && ev.target instanceof Node && !root.current.contains(ev.target)) setOpen(false);
    };
    document.addEventListener("mousedown", onDown, true);
    return () => document.removeEventListener("mousedown", onDown, true);
  }, [open]);
  const close = (): void => {
    setOpen(false);
    button.current?.focus();
  };
  const pick = (sort: T): void => {
    onSort(sort);
    close();
  };
  const onMenuKey = (e: KeyboardEvent<HTMLDivElement>): void => {
    if (e.key === "Escape" || e.key === "Tab") {
      e.preventDefault();
      e.stopPropagation();
      close();
      return;
    }
    const items = Array.from(root.current?.querySelectorAll<HTMLButtonElement>("[role='menuitemradio']") ?? []);
    const at = items.indexOf(document.activeElement as HTMLButtonElement);
    const step = STEP[e.key];
    const next = e.key === "Home" ? 0 : e.key === "End" ? items.length - 1 : step === undefined || at < 0 ? -1 : (at + step + items.length) % items.length;
    if (next < 0) return;
    e.preventDefault();
    items[next]!.focus();
  };
  return (
    <div ref={root} className="relative shrink-0">
      <Tooltip content={`Sort ${what}`}>
        <button
          ref={button}
          type="button"
          aria-haspopup="menu"
          aria-expanded={open}
          aria-label={`Sort ${what} by ${shown?.name ?? ""}`.trim()}
          data-testid="sample-sort"
          onClick={() => setOpen((o) => !o)}
          className="flex h-[30px] items-center gap-1 rounded-sm border border-hair bg-plate px-2 text-xs text-ink-soft hover:text-ink"
        >
          <SortGlyph />
          <span className="grid">
            {choices.map((c) => (
              <span key={c.sort} aria-hidden="true" className={"col-start-1 row-start-1 whitespace-nowrap" + (c.sort === shown?.sort ? "" : " invisible")}>
                {c.name}
              </span>
            ))}
          </span>
        </button>
      </Tooltip>
      {open && (
        <div
          role="menu"
          aria-label={`Sort ${what} by`}
          data-keys-trap
          data-testid="sample-sort-menu"
          onKeyDown={onMenuKey}
          className="absolute right-0 top-full z-20 mt-1 flex w-56 flex-col rounded-sm border border-hair bg-plate p-1 shadow-lg"
        >
          {choices.map((c) => {
            const on = c.sort === current;
            return (
              <button
                key={c.sort}
                type="button"
                role="menuitemradio"
                aria-checked={on}
                tabIndex={on ? 0 : -1}
                data-sort={c.sort}
                onClick={() => pick(c.sort)}
                className={"flex items-start gap-1.5 rounded-sm px-2 py-1.5 text-left text-base text-ink hover:bg-paper-sunk focus:bg-paper-sunk focus:outline-none" + (on ? " bg-paper-sunk" : "")}
              >
                <span aria-hidden="true" className="w-3.5 shrink-0 font-semibold">{on ? "✓" : ""}</span>
                <span className="flex flex-col">
                  <span>{c.name}</span>
                  <span className="text-xs text-ink-faint">{c.hint}</span>
                </span>
              </button>
            );
          })}
        </div>
      )}
    </div>
  );
}
