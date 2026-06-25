import { useRef, useState } from "react";
import type { KeyboardEvent } from "react";
import { Input } from "../ui/Input";
import { Popover } from "../ui/Popover";
import type { ValidatePathResponse } from "../../api";

export interface DirectoryPickerFieldProps {
  value: string;
  onChange: (path: string) => void;
  /** Live autocomplete suggestions for the current `value` (caller fetches via
   *  `suggestPaths`). */
  suggestions: ReadonlyArray<string>;
  /** Live validation probe result, or null when not yet probed. */
  validation: ValidatePathResponse | null;
  /** Accessible label for the input. Default "Data directory" (the primary
   *  new-experiment use); set e.g. "Analysis directory" when reused elsewhere. */
  ariaLabel?: string;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * DirectoryPickerField — path input + suggestion list (via Popover) +
 * Tab/↑↓/↵ completion + a validation line. No full file browser (spec §8.7).
 * Composes `Input`, `Popover`, and design-system token utilities only
 * (no inline appearance); the suggestion list is anchored inside Popover,
 * which provides outside-pointerdown dismiss, Escape close, and focus-return
 * to the trigger — the functional gaps of a hand-rolled `<ul role=listbox>`.
 *
 * Popover trigger: the Input itself (its inner `<input>` receives the
 * `onValueChange` + keyboard handlers via `...rest`; the Popover clone adds
 * `aria-haspopup="dialog"` + `aria-expanded` without conflicting with those).
 * The popover stays open whenever `suggestions.length > 0` — Popover's
 * outside-pointerdown close matches the user expectation of clicking away to
 * dismiss. Escape resets focus to the input (Popover handles it).
 */
export function DirectoryPickerField({
  value,
  onChange,
  suggestions,
  validation,
  ariaLabel = "Data directory",
  className = "",
}: DirectoryPickerFieldProps): JSX.Element {
  const [active, setActive] = useState(0);
  // Ref to the inner <input>. Popover's open-effect focuses its panel by
  // default; for this search-first combobox we redirect that focus back to the
  // input (initialFocusRef) so a click that opens the suggestion popover keeps
  // the caret in the field — otherwise typing and Tab-completion land on the
  // panel, not the input.
  const inputRef = useRef<HTMLInputElement>(null);

  const complete = (path: string): void => {
    onChange(path);
    setActive(0);
  };

  // The dropdown shows each suggestion's basename (last path segment) — the only
  // part that differs from what the user typed — rather than the full absolute
  // path, which is long and repeats the common prefix. Selecting still completes
  // the FULL path into the field (complete(s)). Defensive against a trailing
  // slash; falls back to the whole string if there is no segment.
  const basename = (p: string): string => {
    const trimmed = p.replace(/\/+$/, "");
    const seg = trimmed.slice(trimmed.lastIndexOf("/") + 1);
    return seg === "" ? p : seg;
  };

  const onKeyDown = (e: KeyboardEvent<HTMLInputElement>): void => {
    if (suggestions.length === 0) return;
    if (e.key === "Tab") {
      // Tab completes the TOP suggestion (spec §8.7). preventDefault so focus
      // doesn't leave the field on the completing keystroke.
      e.preventDefault();
      complete(suggestions[0]!);
    } else if (e.key === "ArrowDown") {
      e.preventDefault();
      setActive((i) => Math.min(i + 1, suggestions.length - 1));
    } else if (e.key === "ArrowUp") {
      e.preventDefault();
      setActive((i) => Math.max(i - 1, 0));
    } else if (e.key === "Enter") {
      e.preventDefault();
      complete(suggestions[active] ?? suggestions[0]!);
    }
  };

  // The Input is the Popover trigger. Popover's cloneElement adds
  // aria-haspopup/aria-expanded/onClick to it; onValueChange and onKeyDown
  // pass through Input's ...rest spread onto the inner <input>. Popover is
  // self-managed (no controlled `open`): clicking into the field to type opens
  // it, and the list below renders once suggestions arrive — so in the normal
  // click-then-type flow the matches appear as you type. (A purely programmatic
  // fill that never focuses the field would not auto-open it; not a real-use
  // path. 9d.) Popover's outside-pointerdown + Escape handlers close it.
  const hasSuggestions = suggestions.length > 0;

  return (
    <div className={`flex flex-col gap-2 ${className}`.trim()}>
      <Popover
        trigger={
          <Input
            testId="dirpicker-input"
            value={value}
            onValueChange={onChange}
            onKeyDown={onKeyDown}
            inputRef={inputRef}
            mono
            placeholder="Type or paste your experiment directory…"
            aria-label={ariaLabel}
            className="w-full"
          />
        }
        label="Directory suggestions"
        initialFocusRef={inputRef}
        fullWidth
        className="w-full"
      >
        {hasSuggestions && (
          <ul
            data-testid="dirpicker-suggestions"
            role="listbox"
            aria-label="Directory suggestions"
            className="py-1"
          >
            {suggestions.map((s, i) => (
              <li key={s}>
                <button
                  type="button"
                  role="option"
                  aria-selected={i === active}
                  data-active={i === active ? "true" : undefined}
                  onClick={() => complete(s)}
                  title={s}
                  className={
                    "flex w-full px-3 py-1.5 text-left font-mono text-sm transition-colors " +
                    (i === active ? "text-ink bg-paper-sunk" : "text-ink-soft hover:text-ink hover:bg-paper-sunk")
                  }
                >
                  {basename(s)}
                </button>
              </li>
            ))}
          </ul>
        )}
      </Popover>
      {validation && (
        <div
          data-testid="dirpicker-validation"
          data-ok={validation.ok ? "true" : "false"}
          className={"text-sm " + (validation.ok ? "text-ink-soft" : "text-error")}
        >
          {validation.ok
            ? `Matched ${validation.matched} of ${validation.scanned} files`
            : (validation.message ?? "No exposures found")}
        </div>
      )}
    </div>
  );
}
