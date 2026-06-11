import { useEffect, useId, useRef, useState } from "react";
import type { KeyboardEvent } from "react";
import { Input } from "./Input";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface TagSuggestOption {
  label: string;
  count?: number;
}

export interface TagSuggestProps {
  /** Accessible name for the combobox control. */
  label: string;
  /** "key" vs "value" mode (exposed via data-mode for tests). */
  mode: "key" | "value";
  /** Controlled text value. */
  value: string;
  /** Candidate options with optional usage counts. */
  options: TagSuggestOption[];
  /** Called with the committed string (Enter / option click). */
  onCommit: (v: string) => void;
  /** Called on any text change (update the controlled value). */
  onValueChange?: (v: string) => void;
  /** When true, give the inner <input> an error border + aria-invalid. */
  invalid?: boolean;
  /** aria-describedby id for the input when invalid. */
  "aria-describedby"?: string;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * TagSuggest — a WAI-ARIA combobox for both the key and value fields of the
 * ManageTagsModal row editor. Composes the base Input primitive (appearance
 * lives there) and renders a `role="listbox"` panel styled on the Menu/Popover
 * vocabulary (`.card` Plate-Lift surface). Lives in `src/print/ui/` → design-
 * guard EXEMPT (appearance-authoring zone).
 *
 * ARIA: input is `role="combobox"` + `aria-expanded` + `aria-controls` (the
 * listbox id) + `aria-autocomplete="list"` + `aria-activedescendant` (the
 * focused option id). Each option is `role="option"`. Keyboard: ArrowDown/Up
 * moves the active option; Enter commits the active option (or typed text when
 * nothing is active); Escape collapses.
 *
 * A "create-as-typed" trailing option is shown when the text matches no
 * existing option exactly — committing it yields the typed text verbatim.
 * The did-you-mean near-match nudge is deferred.
 */
export function TagSuggest({
  label,
  mode,
  value,
  options,
  onCommit,
  onValueChange,
  invalid = false,
  "aria-describedby": ariaDescribedBy,
  className = "",
}: TagSuggestProps): JSX.Element {
  const listboxId = useId();
  const optionIdPrefix = useId();

  const [open, setOpen] = useState(false);
  const [activeIdx, setActiveIdx] = useState<number | null>(null);
  const inputRef = useRef<HTMLInputElement>(null);

  // Build displayed option list: filtered existing options + optional create row
  const trimmed = value.trim();
  const exactMatch = options.some(
    (o) => o.label.toLowerCase() === trimmed.toLowerCase(),
  );
  const filteredOptions = trimmed
    ? options.filter((o) =>
        o.label.toLowerCase().includes(trimmed.toLowerCase()),
      )
    : options;
  const showCreate = trimmed !== "" && !exactMatch;

  // Total item count = filtered options + possibly one create row
  const totalItems = filteredOptions.length + (showCreate ? 1 : 0);
  const createIdx = showCreate ? filteredOptions.length : -1;

  // Close and reset active when options list changes meaningfully
  useEffect(() => {
    setActiveIdx(null);
  }, [value]);

  const openList = (): void => {
    if (totalItems > 0) setOpen(true);
  };

  const closeList = (): void => {
    setOpen(false);
    setActiveIdx(null);
  };

  const commitValue = (v: string): void => {
    onCommit(v);
    closeList();
    // Keep focus in the input after commit so the modal row stays navigable
    inputRef.current?.focus();
  };

  const handleInputChange = (v: string): void => {
    onValueChange?.(v);
    // Reopen the list when typing (if there are candidates)
    if (totalItems > 0 || v.trim() !== "") setOpen(true);
    else setOpen(false);
  };

  const handleKeyDown = (e: KeyboardEvent<HTMLInputElement>): void => {
    if (e.key === "ArrowDown") {
      e.preventDefault();
      if (!open) {
        openList();
        setActiveIdx(totalItems > 0 ? 0 : null);
        return;
      }
      if (totalItems === 0) return;
      setActiveIdx((prev) => {
        if (prev === null) return 0;
        return (prev + 1) % totalItems;
      });
    } else if (e.key === "ArrowUp") {
      e.preventDefault();
      if (!open) {
        openList();
        setActiveIdx(totalItems > 0 ? totalItems - 1 : null);
        return;
      }
      if (totalItems === 0) return;
      setActiveIdx((prev) => {
        if (prev === null) return totalItems - 1;
        return (prev - 1 + totalItems) % totalItems;
      });
    } else if (e.key === "Enter") {
      e.preventDefault();
      if (open && activeIdx !== null) {
        if (activeIdx === createIdx) {
          commitValue(trimmed);
        } else {
          const opt = filteredOptions[activeIdx];
          if (opt) commitValue(opt.label);
        }
      } else {
        // No active option: commit the typed text (or do nothing on empty)
        if (trimmed) commitValue(trimmed);
      }
    } else if (e.key === "Escape") {
      if (open) {
        e.preventDefault();
        e.stopPropagation(); // Don't let ModalShell's Escape close the modal
        closeList();
      }
      // If already closed, let Escape propagate (modal handles it)
    }
  };

  const handleBlur = (): void => {
    // Small delay so option clicks can fire before we close
    setTimeout(() => {
      closeList();
    }, 150);
  };

  const activeOptionId =
    open && activeIdx !== null
      ? `${optionIdPrefix}-opt-${activeIdx}`
      : undefined;

  return (
    <div
      data-testid="tag-suggest"
      data-mode={mode}
      className={cx("relative inline-flex", className)}
    >
      <Input
        value={value}
        onValueChange={handleInputChange}
        onKeyDown={handleKeyDown}
        onFocus={openList}
        onBlur={handleBlur}
        invalid={invalid}
        inputSize="sm"
        // Fill the TagSuggest root (which the consumer can stretch via
        // className) so the field tracks its container width instead of sitting
        // at content width. The default `inline-flex` root still collapses to
        // content when the consumer passes no width, preserving every existing
        // consumer's layout.
        className="w-full"
        aria-label={label}
        // Combobox ARIA
        role="combobox"
        aria-expanded={open && totalItems > 0}
        aria-controls={listboxId}
        aria-autocomplete="list"
        {...(activeOptionId ? { "aria-activedescendant": activeOptionId } : {})}
        {...(ariaDescribedBy ? { "aria-describedby": ariaDescribedBy } : {})}
        inputRef={inputRef}
      />
      {open && totalItems > 0 && (
        <div
          id={listboxId}
          role="listbox"
          aria-label={label}
          data-testid="tag-suggest-listbox"
          className="card absolute left-0 top-full mt-0.5 z-30 min-w-[160px] py-1"
        >
          {filteredOptions.map((opt, idx) => {
            const optId = `${optionIdPrefix}-opt-${idx}`;
            const isActive = activeIdx === idx;
            const faint = opt.count !== undefined && opt.count <= 1;
            return (
              <div
                key={opt.label}
                id={optId}
                role="option"
                aria-selected={isActive}
                data-testid="tag-suggest-option"
                data-active={isActive ? "true" : undefined}
                onMouseDown={(e) => {
                  // Prevent input blur from firing before click
                  e.preventDefault();
                  commitValue(opt.label);
                }}
                className={cx(
                  "flex items-center justify-between px-3 py-1 text-sm cursor-pointer",
                  isActive
                    ? "bg-paper-sunk text-ink"
                    : "text-ink-soft hover:bg-paper-sunk hover:text-ink",
                )}
              >
                <span>{opt.label}</span>
                {opt.count !== undefined && (
                  <span
                    data-testid="tag-suggest-count"
                    className={cx(
                      "ml-3 text-xs tabular-nums",
                      faint ? "text-ink-faint" : "text-ink-soft",
                    )}
                  >
                    {opt.count}
                  </span>
                )}
              </div>
            );
          })}
          {showCreate && (
            <div
              id={`${optionIdPrefix}-opt-${createIdx}`}
              role="option"
              aria-selected={activeIdx === createIdx}
              data-testid="tag-suggest-create"
              data-active={activeIdx === createIdx ? "true" : undefined}
              onMouseDown={(e) => {
                e.preventDefault();
                commitValue(trimmed);
              }}
              className={cx(
                "flex items-center px-3 py-1 text-sm cursor-pointer italic",
                activeIdx === createIdx
                  ? "bg-paper-sunk text-ink"
                  : "text-ink-soft hover:bg-paper-sunk hover:text-ink",
              )}
            >
              Create "{trimmed}"
            </div>
          )}
        </div>
      )}
    </div>
  );
}
