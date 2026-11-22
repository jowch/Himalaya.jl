export interface SuggestionItem {
  key: string;
  name: string;
  hint?: string | undefined;
  title?: string | undefined;
  strong?: boolean;
}

export interface SuggestionListProps {
  id: string;
  label: string;
  items: ReadonlyArray<SuggestionItem>;
  active: number;
  onPick: (index: number) => void;
  note?: string | undefined;
  testId?: string;
  className?: string;
}

export const suggestionOptionId = (listId: string, index: number): string => `${listId}-${index}`;

export function SuggestionList({ id, label, items, active, onPick, note, testId, className = "" }: SuggestionListProps): JSX.Element {
  return (
    <div className={`card z-30 p-2 ${className}`.trim()} data-testid={testId}>
      <ul id={id} role="listbox" aria-label={label} className="py-1">
        {items.map((item, i) => (
          <li key={item.key}>
            <button
              type="button"
              tabIndex={-1}
              id={suggestionOptionId(id, i)}
              role="option"
              aria-selected={i === active}
              aria-label={item.hint === undefined ? item.name : `${item.name}, ${item.hint}`}
              data-active={i === active ? "true" : undefined}
              onMouseDown={(e) => e.preventDefault()}
              onClick={() => onPick(i)}
              title={item.title}
              className={
                "flex w-full items-baseline gap-3 px-3 py-1.5 text-left text-sm transition-colors " +
                (i === active ? "text-ink bg-paper-sunk" : "text-ink-soft hover:text-ink hover:bg-paper-sunk")
              }
            >
              <span className={"min-w-0 flex-1 truncate" + (item.strong === true ? " text-ink" : "")}>{item.name}</span>
              {item.hint !== undefined && <span className="shrink-0 text-caption text-ink-soft">{item.hint}</span>}
            </button>
          </li>
        ))}
      </ul>
      {note !== undefined && <p data-testid={testId === undefined ? undefined : `${testId}-note`} className="px-3 py-1 text-caption text-ink-soft">{note}</p>}
    </div>
  );
}
