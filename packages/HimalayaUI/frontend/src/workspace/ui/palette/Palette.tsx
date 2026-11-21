import { useEffect, useId, useMemo, useRef, useState, type KeyboardEvent } from "react";
import { Kbd } from "../Kbd";
import { Overlay } from "../Overlay";
import { pageRows, pickable, rootRows, type PaletteForm, type PaletteGroup, type PaletteItem, type PalettePage, type PaletteRow } from "./model";
import { PaletteFooter } from "./PaletteFooter";
import { PaletteFormView } from "./PaletteFormView";

type Level = { query: string; cursor: number } & ({ page: string; form?: undefined } | { form: PaletteForm; page?: undefined });

interface Props {
  groups: readonly PaletteGroup[];
  pages: Readonly<Record<string, PalettePage>>;
  onClose: () => void;
}

export function Palette({ groups, pages, onClose }: Props): JSX.Element {
  const [stack, setStack] = useState<Level[]>([]);
  const [query, setQuery] = useState("");
  const [cursor, setCursor] = useState(0);
  const listId = useId();
  const list = useRef<HTMLUListElement | null>(null);
  const input = useRef<HTMLInputElement | null>(null);
  const live = stack.filter((l) => l.form !== undefined || pages[l.page] !== undefined);
  const top = live[live.length - 1];
  const form = top?.form;
  const page = top?.page !== undefined ? pages[top.page] : undefined;
  const rows = useMemo<PaletteRow[]>(() => (page ? pageRows(page, query) : rootRows(groups, pages, query)), [page, groups, pages, query]);
  const picks = useMemo(() => pickable(rows), [rows]);
  const at = Math.min(cursor, Math.max(0, picks.length - 1));
  const optionId = (i: number): string => `${listId}-${i}`;

  useEffect(() => {
    const el = list.current?.querySelector<HTMLElement>('[aria-selected="true"]');
    el?.scrollIntoView?.({ block: "nearest" });
  }, [at, rows]);
  useEffect(() => input.current?.focus(), [live.length]);

  const open = (next: { page: string } | { form: PaletteForm }): void => {
    setStack([...live, { ...next, query, cursor: at }]);
    setQuery("");
    setCursor(0);
  };
  const back = (): void => {
    const prev = live[live.length - 1];
    if (!prev) return onClose();
    setStack(live.slice(0, -1));
    setQuery(prev.query);
    setCursor(prev.cursor);
  };
  const choose = (item: PaletteItem | undefined): void => {
    if (!item) return;
    if (item.page && pages[item.page]) return open({ page: item.page });
    if (item.form) return open({ form: item.form() });
    if (!item.run) return;
    onClose();
    item.run();
  };
  const step = (d: number, wrap: boolean): void => {
    const n = picks.length;
    if (n === 0) return;
    setCursor(wrap ? (at + d + n) % n : Math.max(0, Math.min(n - 1, at + d)));
  };
  const onKeyDown = (e: KeyboardEvent<HTMLInputElement>): void => {
    if (e.metaKey || e.ctrlKey || e.altKey) return;
    if (e.key === "ArrowDown" || e.key === "ArrowUp") {
      e.preventDefault();
      step(e.key === "ArrowDown" ? 1 : -1, false);
    } else if (e.key === "Tab") {
      e.preventDefault();
      step(e.shiftKey ? -1 : 1, true);
    } else if (e.key === "Enter") {
      e.preventDefault();
      choose(picks[at]?.item);
    } else if (e.key === "Backspace" && query === "" && page) {
      e.preventDefault();
      back();
    }
  };

  const crumbs = live.slice(0, page?.arg ? -1 : undefined).flatMap((l) => (l.page !== undefined ? [pages[l.page]!.title] : []));
  const overlayKey = (key: string): boolean => key === "Escape" && (back(), true);
  if (form) {
    return (
      <Overlay label={form.title} onClose={onClose} closeKeys={[]} onKey={overlayKey} testId="workspace-palette">
        <PaletteFormView key={live.length} form={form} onBack={back} onDone={onClose} />
      </Overlay>
    );
  }
  const field = (
    <input
      ref={input}
      type="text"
      role="combobox"
      autoComplete="off"
      spellCheck={false}
      aria-expanded="true"
      aria-controls={listId}
      aria-autocomplete="list"
      aria-activedescendant={picks.length > 0 ? optionId(at) : undefined}
      aria-label={page ? (page.arg ?? `Filter ${page.title}`) : "Command or sample"}
      value={query}
      onChange={(e) => {
        setQuery(e.target.value);
        setCursor(0);
      }}
      onKeyDown={onKeyDown}
      placeholder={page ? (page.placeholder ?? "Type to filter") : "Type a command or a sample"}
      className="h-full min-w-0 flex-1 bg-transparent text-lg text-ink placeholder:text-ink-faint focus:outline-none"
    />
  );
  let index = -1;
  return (
    <Overlay label={page ? page.title : "Commands"} onClose={onClose} closeKeys={[]} onKey={overlayKey} testId="workspace-palette">
      <div className="flex h-12 items-center gap-2 overflow-hidden border-b border-hair px-3" data-testid="workspace-palette-bar">
        {crumbs.map((t, i) => (
          <span key={i} className="inline-flex h-6 shrink-0 items-center whitespace-nowrap rounded-sm bg-paper-sunk px-2 text-base text-ink-soft">
            {t} ›
          </span>
        ))}
        {page?.arg ? (
          <div className="flex min-w-0 flex-1 items-center gap-1.5">
            <span className="whitespace-nowrap text-lg font-semibold">{page.title}</span>
            <span className="inline-flex h-7 min-w-0 flex-1 items-center gap-1.5 rounded-sm bg-plate px-2 shadow-[inset_0_0_0_1px_var(--color-accent)]">
              <span className="shrink-0 text-xs text-ink-faint">{page.arg}</span>
              {field}
            </span>
          </div>
        ) : (
          field
        )}
      </div>
      <ul ref={list} id={listId} role="listbox" aria-label={page ? page.title : "Matching commands"} className="max-h-[50vh] overflow-y-auto p-1">
        {rows.map((r, ri) => {
          if (r.kind === "head") {
            return (
              <li key={`h-${ri}`} role="presentation" className="px-2 pb-1 pt-2.5 text-xs text-ink-faint">
                {r.title}
              </li>
            );
          }
          const i = ++index;
          const { item, path } = r;
          return (
            <li
              key={`${path ?? ""}${item.id}`}
              id={optionId(i)}
              role="option"
              aria-selected={i === at}
              aria-label={`${path ? `${path} › ` : ""}${item.label}${item.hint ? `, ${item.hint}` : ""}`}
              data-palette-item={item.id}
              onMouseEnter={() => setCursor(i)}
              onMouseDown={(e) => e.preventDefault()}
              onClick={() => choose(item)}
              className={
                "grid h-8 cursor-pointer grid-cols-[minmax(0,1fr)_auto_minmax(32px,auto)] items-center gap-x-3 rounded-sm px-2 text-base " +
                (i === at ? "bg-row-cursor" : "")
              }
            >
              <span className={"flex min-w-0 items-center gap-2 " + (item.danger ? "text-error" : "")}>
                {item.swatch && <span aria-hidden="true" className={"h-2 w-2 shrink-0 rounded-full " + item.swatch} />}
                <span className="truncate">
                  {path && <span className="text-ink-faint">{path} › </span>}
                  {item.label}
                </span>
              </span>
              <span className="max-w-[260px] truncate text-right text-xs text-ink-faint">{item.hint}</span>
              <span className="flex justify-end text-ink-faint">{item.shortcut ? <Kbd>{item.shortcut}</Kbd> : item.page || item.form ? "›" : null}</span>
            </li>
          );
        })}
        {picks.length === 0 && <li className="px-2 py-2 text-xs text-ink-faint">{page?.empty ?? "No match."}</li>}
      </ul>
      <PaletteFooter keys={[["↵", page ? "pick" : "run"], ["↑ ↓", "move"], ["Esc", page ? "back" : "close"]]} />
    </Overlay>
  );
}
