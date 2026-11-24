import { useEffect, useMemo, useRef, useState } from "react";
import type { KeyboardEvent } from "react";
import { keepPreviousData, useQuery } from "@tanstack/react-query";
import * as api from "../api";
import { queryKeys } from "../queries";
import { FolderBrowser } from "../ui/FolderBrowser";
import type { FolderBrowserRow, FolderCrumb } from "../ui/FolderBrowser";
import { useListHighlight } from "../ui/useListHighlight";
import { useFolderDescriptions } from "../hooks/useFolderDescriptions";

export type FolderPickerMode = "beamtime" | "any";

export interface FolderPickerProps {
  start: string;
  highlight?: string | undefined;
  mode: FolderPickerMode;
  onChoose: (path: string) => void;
  onCancel?: (() => void) | undefined;
  autoFocus?: boolean;
  testId?: string;
  className?: string;
}

const LISTING_FRESH_MS = 10_000;

export const joinFolder = (dir: string, name: string): string => (dir.endsWith("/") ? `${dir}${name}` : `${dir}/${name}`);

export function crumbsOf(path: string): FolderCrumb[] {
  const parts = path.split("/").filter((p) => p !== "");
  return [{ label: "/", path: "/" }, ...parts.map((label, i) => ({ label, path: `/${parts.slice(0, i + 1).join("/")}` }))];
}

const plural = (n: number, one: string): string => `${n} ${one}${n === 1 ? "" : "s"}`;

export function describeHint(d: api.FolderDescription | undefined): string {
  if (d === undefined) return "Checking…";
  const kind = d.kind === "beamtime" ? `Beamtime · ${plural(d.setups ?? 0, "setup")}` : d.kind === "setup" ? "Setup" : "";
  const added = d.added === null ? "" : `added as ${d.added.name}`;
  return [kind, added].filter((s) => s !== "").join(" · ");
}

const isTypedChar = (e: KeyboardEvent): boolean =>
  e.key.length === 1 && e.key !== "/" && e.key !== "?" && e.key !== " " && !e.metaKey && !e.ctrlKey && !e.altKey && !e.nativeEvent.isComposing;

export const folderAbove = (path: string): string => {
  const t = path.replace(/\/+$/, "");
  if (!t.startsWith("/") && !t.startsWith("~")) return "";
  const cut = t.lastIndexOf("/");
  return cut <= 0 ? (t.startsWith("/") ? "/" : "") : t.slice(0, cut);
};

export function FolderPicker({ start, highlight, mode, onChoose, onCancel, autoFocus = false, testId = "folder-browser", className }: FolderPickerProps): JSX.Element {
  const [dir, setDir] = useState(start);
  const [filter, setFilter] = useState("");
  const [editing, setEditing] = useState(false);
  const [draft, setDraft] = useState("");
  const [draftError, setDraftError] = useState<string | undefined>(undefined);
  const [landOn, setLandOn] = useState<string | undefined>(highlight);
  const listRef = useRef<HTMLDivElement>(null);

  const listing = useQuery({
    queryKey: queryKeys.folderListing(dir),
    queryFn: ({ signal }) => api.listFolder(dir, { signal }),
    staleTime: LISTING_FRESH_MS,
    retry: false,
    placeholderData: keepPreviousData,
  });
  const here = listing.data;
  const lost = listing.isError && here === undefined && dir !== "" && dir !== "/";
  useEffect(() => {
    if (lost) setDir(folderAbove(dir));
  }, [lost, dir]);
  const shown = useMemo(() => {
    if (here === undefined) return [];
    const needle = filter.toLowerCase();
    return here.folders.filter((n) => n.toLowerCase().includes(needle)).map((n) => joinFolder(here.path, n));
  }, [here, filter]);
  const describing = mode === "beamtime" && !listing.isPlaceholderData && here?.inside_setup !== true;
  const { found, describeNow } = useFolderDescriptions(shown, describing);

  const choose = async (path: string): Promise<void> => {
    if (mode === "any") { onChoose(path); return; }
    if (here?.inside_setup === true) { open(path); return; }
    const d = found.get(path) ?? await describeNow(path).catch(() => undefined);
    if (d !== undefined && d.kind !== null) onChoose(d.path);
    else open(path);
  };
  const open = (path: string): void => {
    setDir(path);
    setFilter("");
  };
  const up = (): void => {
    if (here?.parent == null) return;
    setLandOn(here.path);
    open(here.parent);
  };

  const hl = useListHighlight({
    count: shown.length,
    keys: shown,
    resetKey: `${here?.path ?? ""}\n${filter}`,
    onCommit: (i) => { const p = shown[i]; if (p !== undefined) void choose(p); },
    onEscape: () => (filter !== "" ? setFilter("") : onCancel !== undefined ? onCancel() : up()),
  });
  const { setAt } = hl;

  useEffect(() => {
    if (landOn === undefined || here === undefined || listing.isPlaceholderData) return;
    const i = shown.indexOf(landOn);
    if (i >= 0) setAt(i);
    setLandOn(undefined);
  }, [landOn, here, shown, listing.isPlaceholderData, setAt]);

  useEffect(() => {
    if (autoFocus) listRef.current?.focus();
  }, [autoFocus]);

  const startEditing = (): void => {
    setDraft(here === undefined ? "" : joinFolder(here.path, ""));
    setDraftError(undefined);
    setEditing(true);
  };
  const stopEditing = (): void => {
    setEditing(false);
    setDraftError(undefined);
    requestAnimationFrame(() => listRef.current?.focus());
  };
  const submitDraft = async (): Promise<void> => {
    const typed = draft.trim();
    if (typed === "") return;
    try {
      const target = await api.listFolder(typed);
      stopEditing();
      if (mode === "any") { onChoose(target.path); return; }
      const d = await describeNow(target.path).catch(() => undefined);
      if (d !== undefined && d.kind !== null) onChoose(d.path);
      else open(target.path);
    } catch (err) {
      setDraftError(err instanceof api.ApiError ? sentence(err.message) : "Could not read that folder.");
    }
  };

  const onListKeyDown = (e: KeyboardEvent<HTMLDivElement>): void => {
    if (e.metaKey || e.ctrlKey || e.altKey) return;
    if (e.key === "ArrowRight") {
      const p = shown[hl.at];
      if (p !== undefined) { e.preventDefault(); open(p); }
      return;
    }
    if (e.key === "ArrowLeft") { e.preventDefault(); up(); return; }
    if (e.key === "/") { e.preventDefault(); startEditing(); return; }
    if (e.key === "Backspace") {
      if (filter !== "") { e.preventDefault(); setFilter((f) => f.slice(0, -1)); }
      return;
    }
    if (hl.onKeyDown(e)) return;
    if (isTypedChar(e)) { e.preventDefault(); setFilter((f) => f + e.key); }
  };

  const onDraftKeyDown = (e: KeyboardEvent<HTMLInputElement>): void => {
    if (e.nativeEvent.isComposing) return;
    if (e.key === "Enter") { e.preventDefault(); void submitDraft(); }
    if (e.key === "Escape") { e.preventDefault(); e.stopPropagation(); stopEditing(); }
  };

  const rows: FolderBrowserRow[] = shown.map((path) => {
    const d = found.get(path);
    const name = path.slice(path.lastIndexOf("/") + 1);
    return describing
      ? { path, name, hint: describeHint(d), waiting: d === undefined, strong: d !== undefined && d.kind !== null }
      : { path, name, hint: "", waiting: false, strong: false };
  });

  const status = lost
    ? "Opening folder…"
    : listing.isError
    ? listing.error instanceof api.ApiError ? sentence(listing.error.message) : "Could not read this folder."
    : here === undefined || listing.isPlaceholderData
      ? "Opening folder…"
      : here.folders.length === 0
        ? "No folders here."
        : shown.length === 0
          ? `No folder here has “${filter}” in its name.`
          : undefined;

  return (
    <FolderBrowser
      testId={testId}
      className={className}
      crumbs={here === undefined ? [] : crumbsOf(here.path)}
      rows={rows}
      active={hl.at}
      filter={filter}
      status={status}
      editing={editing}
      draft={draft}
      draftError={draftError}
      chooseHint={mode === "any" ? "to use this folder" : "to add a beamtime or setup, or open a folder"}
      escCloses={onCancel !== undefined}
      listRef={listRef}
      onCrumb={(p) => { open(p); listRef.current?.focus(); }}
      onPick={(i) => { const p = shown[i]; if (p !== undefined) void choose(p); }}
      onOpen={(i) => { const p = shown[i]; if (p !== undefined) open(p); listRef.current?.focus(); }}
      onHover={setAt}
      onListKeyDown={onListKeyDown}
      onStartEditing={startEditing}
      onDraftChange={(v) => { setDraft(v); setDraftError(undefined); }}
      onDraftKeyDown={onDraftKeyDown}
    />
  );
}

function sentence(message: string): string {
  const s = message.trim();
  if (s === "") return s;
  const capital = s[0]!.toUpperCase() + s.slice(1);
  return /[.!?]$/.test(capital) ? capital : `${capital}.`;
}
