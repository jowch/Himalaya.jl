import { useEffect, useId, useMemo, useRef, useState } from "react";
import { Kbd } from "./Kbd";
import { Overlay } from "./Overlay";

export interface Command {
  id: string;
  label: string;
  hint?: string;
  shortcut?: string;
  run: () => void;
}

export function filterCommands(commands: readonly Command[], query: string): Command[] {
  const words = query.toLowerCase().split(/\s+/).filter(Boolean);
  return commands.filter((c) => words.every((w) => c.label.toLowerCase().includes(w)));
}

export function CommandPalette({ commands, onClose }: { commands: readonly Command[]; onClose: () => void }): JSX.Element {
  const [query, setQuery] = useState("");
  const [cursor, setCursor] = useState(0);
  const listId = useId();
  const list = useRef<HTMLUListElement | null>(null);
  const shown = useMemo(() => filterCommands(commands, query).slice(0, 40), [commands, query]);
  const at = Math.min(cursor, Math.max(0, shown.length - 1));
  const optionId = (i: number): string => `${listId}-${i}`;
  useEffect(() => {
    const el = list.current?.querySelector<HTMLElement>('[aria-selected="true"]');
    el?.scrollIntoView?.({ block: "nearest" });
  }, [at, shown]);
  const run = (c: Command | undefined): void => {
    if (!c) return;
    onClose();
    c.run();
  };
  return (
    <Overlay label="Commands" onClose={onClose} testId="workspace-palette">
      <input
        type="text"
        role="combobox"
        aria-expanded="true"
        aria-controls={listId}
        aria-autocomplete="list"
        aria-activedescendant={shown.length > 0 ? optionId(at) : undefined}
        value={query}
        onChange={(e) => {
          setQuery(e.target.value);
          setCursor(0);
        }}
        onKeyDown={(e) => {
          if (e.key === "ArrowDown") {
            e.preventDefault();
            setCursor(Math.min(shown.length - 1, at + 1));
          } else if (e.key === "ArrowUp") {
            e.preventDefault();
            setCursor(Math.max(0, at - 1));
          } else if (e.key === "Enter") {
            e.preventDefault();
            run(shown[at]);
          }
        }}
        placeholder="Type a command or a sample name"
        aria-label="Command or sample"
        className="w-full border-b border-hair bg-plate px-4 py-3 text-lg placeholder:text-ink-faint focus:outline-none"
      />
      <ul ref={list} id={listId} className="max-h-[50vh] overflow-y-auto p-1" role="listbox" aria-label="Matching commands">
        {shown.map((c, i) => (
          <li
            key={c.id}
            id={optionId(i)}
            role="option"
            aria-selected={i === at}
            onMouseEnter={() => setCursor(i)}
            onMouseDown={(e) => e.preventDefault()}
            onClick={() => run(c)}
            className={"flex w-full cursor-pointer justify-between gap-3 rounded-sm px-3 py-1.5 text-left text-body " + (i === at ? "bg-row-cursor" : "")}
          >
            <span>{c.label}</span>
            {c.shortcut ? <Kbd>{c.shortcut}</Kbd> : c.hint && <span className="text-caption">{c.hint}</span>}
          </li>
        ))}
        {shown.length === 0 && <li className="px-3 py-2 text-caption">No match.</li>}
      </ul>
    </Overlay>
  );
}
