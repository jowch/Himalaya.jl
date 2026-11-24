import type { KeyboardEvent, Ref } from "react";
import { cx } from "../lib/cx";
import { Input } from "./Input";
import { Kbd } from "./Kbd";
import { SELECTED_ROW } from "./selectedRow";

export interface FolderCrumb { label: string; path: string }

export interface FolderBrowserRow {
  path: string;
  name: string;
  hint: string;
  waiting: boolean;
  strong: boolean;
}

export interface FolderBrowserProps {
  crumbs: ReadonlyArray<FolderCrumb>;
  rows: ReadonlyArray<FolderBrowserRow>;
  active: number;
  filter: string;
  status: string | undefined;
  editing: boolean;
  draft: string;
  draftError: string | undefined;
  chooseHint: string;
  listRef: Ref<HTMLDivElement>;
  testId?: string;
  className?: string | undefined;
  onCrumb: (path: string) => void;
  onPick: (index: number) => void;
  onOpen: (index: number) => void;
  onHover: (index: number) => void;
  onListKeyDown: (e: KeyboardEvent<HTMLDivElement>) => void;
  onStartEditing: () => void;
  onDraftChange: (v: string) => void;
  onDraftKeyDown: (e: KeyboardEvent<HTMLInputElement>) => void;
}

export const folderRowId = (testId: string, index: number): string => `${testId}-row-${index}`;

export function FolderBrowser({
  crumbs, rows, active, filter, status, editing, draft, draftError, chooseHint, listRef,
  testId = "folder-browser", className = "",
  onCrumb, onPick, onOpen, onHover, onListKeyDown, onStartEditing, onDraftChange, onDraftKeyDown,
}: FolderBrowserProps): JSX.Element {
  const here = crumbs[crumbs.length - 1];
  return (
    <div data-testid={testId} className={cx("flex flex-col gap-2", className)}>
      {editing ? (
        <div className="flex flex-col gap-1">
          <Input
            testId={`${testId}-path`}
            value={draft}
            onValueChange={onDraftChange}
            onKeyDown={onDraftKeyDown}
            aria-label="Folder path"
            placeholder="A full path, starting with / or ~"
            invalid={draftError !== undefined}
            autoFocus
            spellCheck={false}
            autoComplete="off"
            className="w-full"
          />
          {draftError !== undefined && <p role="alert" className="text-meta text-error">{draftError}</p>}
        </div>
      ) : (
        <div data-testid={`${testId}-location`} className="flex flex-wrap items-center gap-1 rounded-sm border border-hair-strong bg-plate px-2.5 py-1.5">
          <nav aria-label="Location" className="flex min-w-0 flex-wrap items-center gap-1">
            {crumbs.map((c, i) => (
              <span key={c.path} className="flex items-center gap-1">
                {i > 0 && <span aria-hidden className="text-ink-faint">›</span>}
                <button
                  type="button"
                  tabIndex={-1}
                  onClick={() => onCrumb(c.path)}
                  aria-current={c === here ? "location" : undefined}
                  className={cx(
                    "rounded-sm px-1.5 transition-colors hover:bg-paper-sunk",
                    c === here ? "font-semibold text-ink" : "text-ink-soft hover:text-ink",
                  )}
                >
                  {c.label}
                </button>
              </span>
            ))}
            {filter !== "" && (
              <span data-testid={`${testId}-filter`} className="ml-2 text-meta text-ink-soft">Showing names with “{filter}”</span>
            )}
          </nav>
          <button
            type="button"
            tabIndex={-1}
            onClick={onStartEditing}
            className="ml-auto flex items-center gap-1.5 text-meta text-ink-soft transition-colors hover:text-ink"
          >
            <Kbd>/</Kbd> Type or paste a path
          </button>
        </div>
      )}
      <div
        ref={listRef}
        role="listbox"
        tabIndex={0}
        aria-label={here === undefined ? "Folders" : `Folders in ${here.path}`}
        aria-activedescendant={rows.length > 0 ? folderRowId(testId, active) : undefined}
        onKeyDown={onListKeyDown}
        data-testid={`${testId}-list`}
        className="max-h-80 overflow-y-auto rounded-sm border border-hair focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent"
      >
        {rows.map((r, i) => (
          <div
            key={r.path}
            id={folderRowId(testId, i)}
            role="option"
            aria-selected={i === active}
            aria-label={r.hint === "" ? r.name : `${r.name}, ${r.hint}`}
            data-testid={`${testId}-row`}
            data-path={r.path}
            onMouseMove={() => onHover(i)}
            onClick={() => onPick(i)}
            className={cx(
              "flex cursor-pointer items-center gap-3 border-t border-hair px-3 py-1.5 first:border-t-0",
              i === active ? SELECTED_ROW : "hover:bg-paper-sunk",
            )}
          >
            <FolderIcon />
            <span className={cx("min-w-0 flex-1 truncate", r.strong ? "font-semibold text-ink" : "text-ink")}>{r.name}</span>
            {r.hint !== "" && (
              <span className={cx("shrink-0 text-meta", r.waiting ? "text-ink-faint" : "text-ink-soft")}>{r.hint}</span>
            )}
            <button
              type="button"
              tabIndex={-1}
              aria-label={`Open ${r.name}`}
              onClick={(e) => { e.stopPropagation(); onOpen(i); }}
              className={cx("shrink-0 rounded-sm px-1 text-ink-faint transition-colors hover:text-ink", i === active ? "visible" : "invisible")}
            >
              ›
            </button>
          </div>
        ))}
      </div>
      {status !== undefined && <p data-testid={`${testId}-status`} role="status" className="text-meta text-ink-soft">{status}</p>}
      <div className="flex flex-wrap gap-x-4 gap-y-1 text-meta text-ink-soft">
        <span className="flex items-center gap-1"><Kbd>↑</Kbd><Kbd>↓</Kbd> or <Kbd>Tab</Kbd> to move</span>
        <span className="flex items-center gap-1"><Kbd>↵</Kbd> {chooseHint}</span>
        <span className="flex items-center gap-1"><Kbd>→</Kbd> to open a folder</span>
        <span className="flex items-center gap-1"><Kbd>←</Kbd> or <Kbd>Esc</Kbd> to go up one folder</span>
      </div>
    </div>
  );
}

function FolderIcon(): JSX.Element {
  return (
    <svg aria-hidden viewBox="0 0 16 13" className="h-3 w-4 shrink-0 text-ink-faint" fill="none" stroke="currentColor" strokeWidth="1.4">
      <path d="M1 2.2c0-.6.4-1 1-1h3.6l1.4 1.6H14c.6 0 1 .4 1 1v7.4c0 .6-.4 1-1 1H2c-.6 0-1-.4-1-1z" />
    </svg>
  );
}
