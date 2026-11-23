import { useEffect, useId, useMemo, useRef, useState, type KeyboardEvent, type MouseEvent } from "react";
import { Kbd } from "../Kbd";
import { Overlay } from "../Overlay";
import { useListHighlight } from "../useListHighlight";
import { pageRows, pickable, rootRows, sections, type PaletteForm, type PaletteGroup, type PaletteItem, type PalettePage, type PaletteRow } from "./model";
import { PaletteFooter } from "./PaletteFooter";
import { PaletteFormView } from "./PaletteFormView";

type Level = { query: string; cursor: number } & ({ page: string; form?: undefined } | { form: PaletteForm; page?: undefined });

interface Props {
  groups: readonly PaletteGroup[];
  pages: Readonly<Record<string, PalettePage>>;
  onClose: () => void;
  onType?: (query: string) => void;
  start?: PaletteForm | undefined;
}

export function keepFocusInField(e: MouseEvent): void {
  if (!(e.target instanceof HTMLInputElement)) e.preventDefault();
}

export function Palette({ groups, pages, onClose, onType, start }: Props): JSX.Element {
  const [stack, setStack] = useState<Level[]>(() => (start ? [{ form: start, query: "", cursor: 0 }] : []));
  const [query, setQuery] = useState("");
  const [ticked, setTicked] = useState<readonly string[]>([]);
  const listId = useId();
  const list = useRef<HTMLDivElement | null>(null);
  const input = useRef<HTMLInputElement | null>(null);
  const heldBack = useRef(false);
  const top = stack[stack.length - 1];
  const form = top?.form;
  const page = top?.page !== undefined ? pages[top.page] : undefined;
  const rows = useMemo<PaletteRow[]>(() => (page ? pageRows(page, query) : rootRows(groups, pages, query)), [page, groups, pages, query]);
  const picks = useMemo(() => pickable(rows), [rows]);
  const optionId = (i: number): string => `${listId}-${i}`;

  const back = (): void => {
    const prev = stack[stack.length - 1];
    if (!prev || (start && stack.length === 1)) return onClose();
    heldBack.current = true;
    if (prev.page !== undefined) setTicked([]);
    setStack(stack.slice(0, -1));
    setQuery(prev.query);
    highlight.setAt(prev.cursor);
  };
  const open = (next: { page: string } | { form: PaletteForm }): void => {
    if ("page" in next) setTicked([]);
    setStack([...stack, { ...next, query, cursor: highlight.at }]);
    setQuery("");
    highlight.setAt(0);
  };
  const multi = page?.multi;
  const tick = (item: PaletteItem): void => setTicked((t) => (t.includes(item.id) ? t.filter((x) => x !== item.id) : [...t, item.id]));
  const commitTicked = (): void => {
    if (!multi) return;
    const next = multi.commit(ticked);
    setTicked([]);
    if (next) return open({ form: next });
    onClose();
  };
  const choose = (item: PaletteItem | undefined): void => {
    if (!item) return;
    if (multi && ticked.length > 0) return commitTicked();
    if (item.page && pages[item.page]) return open({ page: item.page });
    if (item.form) return open({ form: item.form() });
    if (!item.run) return;
    onClose();
    item.run();
  };
  const highlight = useListHighlight({ count: picks.length, onCommit: (i) => choose(picks[i]?.item), onEscape: back, enabled: !form });
  const at = highlight.at;

  useEffect(() => {
    const gone = stack.findIndex((l) => l.page !== undefined && pages[l.page] === undefined);
    if (gone < 0) return;
    setStack(stack.slice(0, gone));
    setQuery(stack[gone]!.query);
    highlight.setAt(stack[gone]!.cursor);
  }, [stack, pages, highlight]);
  useEffect(() => {
    const el = list.current?.querySelector<HTMLElement>('[aria-selected="true"]');
    el?.scrollIntoView?.({ block: "nearest" });
  }, [at, rows]);
  useEffect(() => input.current?.focus(), [stack.length]);

  const onKeyDown = (e: KeyboardEvent<HTMLInputElement>): void => {
    if (e.key === "Backspace" && e.repeat && heldBack.current) return e.preventDefault();
    if (multi && e.key === "Enter" && e.shiftKey && !e.metaKey && !e.ctrlKey && !e.altKey) {
      e.preventDefault();
      const item = picks[highlight.at]?.item;
      if (item && !e.repeat) tick(item);
      return;
    }
    if (highlight.onKeyDown(e)) return;
    if (e.key === "Backspace" && query === "" && page && !e.metaKey && !e.ctrlKey && !e.altKey) {
      e.preventDefault();
      if (!e.repeat) back();
    }
  };

  const crumbs = stack.slice(0, page?.arg ? -1 : undefined).flatMap((l) => (l.page !== undefined && pages[l.page] ? [pages[l.page]!.title] : []));
  if (form) {
    return (
      <Overlay label={form.title} onClose={onClose} closeKeys={[]} testId="workspace-palette">
        <div onMouseDown={keepFocusInField}>
          <PaletteFormView key={stack.length} form={form} onBack={back} onDone={onClose} />
        </div>
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
        onType?.(e.target.value);
        highlight.setAt(0);
      }}
      onKeyDown={onKeyDown}
      onKeyUp={() => {
        heldBack.current = false;
      }}
      placeholder={page ? (page.placeholder ?? "Type to filter") : "Type a command or a sample"}
      className="h-full min-w-0 flex-1 bg-transparent text-lg text-ink placeholder:text-ink-faint focus:outline-none"
    />
  );
  return (
    <Overlay label={page ? page.title : "Commands"} onClose={onClose} closeKeys={[]} testId="workspace-palette">
      <div onMouseDown={keepFocusInField}>
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
        <div ref={list} className="max-h-[50vh] overflow-y-auto p-1">
          <div id={listId} role="listbox" aria-label={page ? page.title : "Matching commands"} aria-multiselectable={multi ? true : undefined}>
            {sections(rows).map((sec, si) => (
              <div key={si} role="group" aria-labelledby={sec.title ? `${listId}-g${si}` : undefined} aria-label={sec.title ? undefined : page?.title}>
                {sec.title && (
                  <div id={`${listId}-g${si}`} role="presentation" className="px-2 pb-1 pt-2.5 text-xs text-ink-faint">
                    {sec.title}
                  </div>
                )}
                {sec.rows.map(({ item, path, index: i }) => (
                  <div
                    key={`${path ?? ""}${item.id}`}
                    id={optionId(i)}
                    role="option"
                    aria-selected={i === at}
                    aria-checked={multi ? ticked.includes(item.id) : undefined}
                    aria-label={`${path ? `${path} › ` : ""}${item.label}${item.value !== undefined ? `, ${item.value || "None"}` : ""}${item.hint ? `, ${item.hint}` : ""}`}
                    data-palette-item={item.id}
                    onMouseEnter={() => highlight.setAt(i)}
                    onMouseDown={(e) => e.preventDefault()}
                    onClick={() => (multi && ticked.length > 0 ? tick(item) : choose(item))}
                    className={
                      "grid h-8 cursor-pointer items-center gap-x-3 rounded-sm px-2 text-base " +
                      (item.value !== undefined ? "grid-cols-[140px_minmax(0,1fr)_auto_minmax(32px,auto)] " : "grid-cols-[minmax(0,1fr)_auto_minmax(32px,auto)] ") +
                      (i === at ? "bg-row-cursor" : "")
                    }
                  >
                    {item.value !== undefined ? (
                      <>
                        <span className="truncate text-ink-soft">{item.label}</span>
                        <span className={"truncate " + (item.value ? "" : "text-ink-faint")}>{item.value || "None"}</span>
                      </>
                    ) : (
                      <span className={"flex min-w-0 items-center gap-2 " + (item.danger ? "text-error" : "")}>
                        {multi && (
                          <span
                            data-tick
                            aria-hidden="true"
                            onClick={(e) => {
                              e.stopPropagation();
                              tick(item);
                            }}
                            className={"inline-flex h-3.5 w-3.5 shrink-0 items-center justify-center rounded-[3px] border text-[10px] leading-none " + (ticked.includes(item.id) ? "border-accent bg-accent text-paper" : "border-hair-strong")}
                          >
                            {ticked.includes(item.id) ? "✓" : null}
                          </span>
                        )}
                        {item.swatch && <span aria-hidden="true" className={"h-2 w-2 shrink-0 rounded-full " + item.swatch} />}
                        <span className="truncate">
                          {path && <span className="text-ink-faint">{path} › </span>}
                          {item.label}
                        </span>
                      </span>
                    )}
                    <span className="max-w-[260px] truncate text-right text-xs text-ink-faint">{item.hint}</span>
                    <span className="flex justify-end text-ink-faint">{item.shortcut ? <Kbd>{item.shortcut}</Kbd> : item.page || item.form ? "›" : null}</span>
                  </div>
                ))}
              </div>
            ))}
          </div>
          {picks.length === 0 && <p className="px-2 py-2 text-xs text-ink-faint">{page?.empty ?? "No match."}</p>}
        </div>
        <PaletteFooter
          keys={[
            ["↵", multi && ticked.length > 0 ? multi.verb(ticked.length) : page ? "pick" : "run"],
            ...(multi ? [["⇧ ↵", "tick"] as const] : []),
            ["Tab", "move"],
            ["Esc", page ? "back" : "close"],
          ]}
        />
      </div>
    </Overlay>
  );
}
