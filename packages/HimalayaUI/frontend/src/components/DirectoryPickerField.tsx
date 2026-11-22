import { useEffect, useId, useRef, useState } from "react";
import type { KeyboardEvent } from "react";
import type { FolderSuggestion } from "../api";
import { Input } from "../ui/Input";
import { SuggestionList, suggestionOptionId } from "../ui/SuggestionList";

export interface DirectoryPickerFieldProps {
  value: string;
  onChange: (path: string) => void;
  onSubmit?: (() => void) | undefined;
  folders: ReadonlyArray<FolderSuggestion>;
  more?: boolean;
  ariaLabel?: string;
  placeholder?: string;
  inputSize?: "sm" | "md";
  testId?: string;
  autoFocus?: boolean;
  className?: string;
}

const trimSlash = (p: string): string => p.replace(/\/+$/, "");

const basename = (p: string): string => {
  const trimmed = trimSlash(p);
  const seg = trimmed.slice(trimmed.lastIndexOf("/") + 1);
  return seg === "" ? p : seg;
};

export function folderHint(f: FolderSuggestion): string | undefined {
  const parts: string[] = [];
  if (f.kind === "beamtime") parts.push(f.setups === 1 ? "beamtime, 1 setup" : `beamtime, ${f.setups ?? 0} setups`);
  if (f.kind === "setup") parts.push("setup");
  if (f.added !== null) parts.push(`added as ${f.added.name}`);
  return parts.length > 0 ? parts.join(" · ") : undefined;
}

export function DirectoryPickerField({
  value,
  onChange,
  onSubmit,
  folders,
  more = false,
  ariaLabel = "Data directory",
  placeholder = "Type or paste a beamtime folder, starting with / or ~",
  inputSize = "md",
  testId = "dirpicker",
  autoFocus = false,
  className = "",
}: DirectoryPickerFieldProps): JSX.Element {
  const [active, setActive] = useState(-1);
  const [focused, setFocused] = useState(autoFocus);
  const [dismissed, setDismissed] = useState(false);
  const listId = useId();
  const root = useRef<HTMLDivElement>(null);
  const foldersKey = folders.map((f) => f.path).join("\n");
  useEffect(() => setActive(-1), [foldersKey]);

  const onlyEcho = folders.length === 1 && trimSlash(folders[0]!.path).toLowerCase() === trimSlash(value).toLowerCase() && !value.endsWith("/");
  const showing = focused && !dismissed && folders.length > 0 && !onlyEcho;
  const highlighted = showing && active >= 0 ? folders[active] : undefined;

  const change = (v: string): void => {
    setDismissed(false);
    onChange(v);
  };

  const complete = (path: string): void => {
    change(path.endsWith("/") ? path : `${path}/`);
    setActive(-1);
  };

  const step = (dir: 1 | -1): void => {
    const n = folders.length;
    setActive((i) => (i < 0 ? (dir === 1 ? 0 : n - 1) : (i + dir + n) % n));
  };

  const onKeyDown = (e: KeyboardEvent<HTMLInputElement>): void => {
    if (e.nativeEvent.isComposing || e.metaKey || e.ctrlKey || e.altKey) return;
    if (e.key === "Escape") {
      if (!showing) return;
      e.preventDefault();
      e.stopPropagation();
      if (active >= 0) setActive(-1);
      else setDismissed(true);
      return;
    }
    if (e.key === "Enter") {
      if (e.shiftKey) return;
      e.preventDefault();
      if (highlighted !== undefined) complete(highlighted.path);
      else onSubmit?.();
      return;
    }
    const down = e.key === "ArrowDown" || (e.key === "Tab" && !e.shiftKey);
    const up = e.key === "ArrowUp" || (e.key === "Tab" && e.shiftKey);
    if (!down && !up) return;
    if (!showing) {
      if (e.key === "ArrowDown" && folders.length > 0) {
        e.preventDefault();
        setDismissed(false);
        setActive(0);
      }
      return;
    }
    if (e.key === "Tab" && e.shiftKey && active < 0) return;
    e.preventDefault();
    step(down ? 1 : -1);
  };

  return (
    <div
      ref={root}
      className={`relative ${className}`.trim()}
      onFocus={() => setFocused(true)}
      onBlur={(e) => {
        if (!(e.relatedTarget instanceof Node && root.current?.contains(e.relatedTarget))) setFocused(false);
      }}
    >
      <Input
        testId={`${testId}-input`}
        inputSize={inputSize}
        value={value}
        onValueChange={change}
        onKeyDown={onKeyDown}
        placeholder={placeholder}
        role="combobox"
        aria-label={ariaLabel}
        aria-autocomplete="list"
        aria-expanded={showing}
        aria-controls={showing ? listId : undefined}
        aria-activedescendant={highlighted === undefined ? undefined : suggestionOptionId(listId, active)}
        autoFocus={autoFocus}
        spellCheck={false}
        autoComplete="off"
        className="w-full"
      />
      {showing && (
        <SuggestionList
          id={listId}
          label="Folders"
          testId={`${testId}-suggestions`}
          className="absolute left-0 right-0 top-full mt-1"
          items={folders.map((f) => ({ key: f.path, name: basename(f.path), hint: folderHint(f), title: f.path, strong: f.kind !== null }))}
          active={active}
          onPick={(i) => complete(folders[i]!.path)}
          note={more ? "More folders match. Keep typing to narrow the list." : undefined}
        />
      )}
    </div>
  );
}
