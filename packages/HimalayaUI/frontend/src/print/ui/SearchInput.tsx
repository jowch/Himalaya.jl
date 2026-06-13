import { useRef } from "react";
import { Input } from "./Input";

interface SearchInputProps {
  value: string;
  onChange: (v: string) => void;
  placeholder?: string;
  /** Accessible name for the field. A placeholder is NOT a reliable label
   *  (WCAG 3.3.2 — it vanishes on input and is only a last-resort accname),
   *  so callers pass an explicit label that wins over the placeholder. */
  ariaLabel?: string;
  className?: string;
}

/** The leading magnifier glyph. Stroke is sourced from a token
 *  (`var(--color-ink-faint)`), never hex; the svg is decorative so it carries
 *  `aria-hidden`. */
function MagnifierSvg(): JSX.Element {
  return (
    <svg
      width={14}
      height={14}
      viewBox="0 0 14 14"
      fill="none"
      aria-hidden="true"
      className="flex-shrink-0"
    >
      <circle cx="6" cy="6" r="4.3" stroke="var(--color-ink-faint)" strokeWidth={1.5} />
      <line
        x1="9.2"
        y1="9.2"
        x2="12.5"
        y2="12.5"
        stroke="var(--color-ink-faint)"
        strokeWidth={1.5}
        strokeLinecap="round"
      />
    </svg>
  );
}

/** Standard product search field: composes the base {@link Input} with a leading
 *  magnifier glyph and, while non-empty, a trailing one-click clear (×).
 *  Appearance (plate well, hairline, 5px radius, the rationed
 *  `focus-within:border-accent` ring) lives in `Input`; SearchInput contributes
 *  only the adornments and keeps its own `data-testid="search-input"` contract
 *  via Input's `testId` override.
 *
 *  The clear × is the house ghost-icon idiom (Chip's remove precedent): a real
 *  in-flow `h-6 w-6` (24×24 CSS px, WCAG 2.5.8 / LO-TAGTARGET) border box whose
 *  negative margins collapse the LAYOUT footprint so the field never inflates
 *  when the button appears. Clearing keeps focus in the input (the user is
 *  mid-search, not done searching).
 *
 *  C — focus is Input's `focus-within:border-accent` ring; hover N/A (a field,
 *  not a button). E — input text is sans/base (prose entry), correct. F —
 *  magnifier stroke sourced from a token (`var(--color-ink-faint)`), never
 *  hex. */
export function SearchInput({
  value,
  onChange,
  placeholder,
  ariaLabel,
  className,
}: SearchInputProps): JSX.Element {
  const inputRef = useRef<HTMLInputElement>(null);
  const clearButton = (
    <button
      type="button"
      aria-label="Clear search"
      onClick={() => {
        onChange("");
        inputRef.current?.focus();
      }}
      className="inline-flex h-6 w-6 shrink-0 items-center justify-center -my-1 -mr-1 rounded-sm text-ink-faint hover:text-ink focus-visible:outline focus-visible:outline-1 focus-visible:outline-offset-0 focus-visible:outline-accent"
    >
      ×
    </button>
  );
  return (
    <Input
      testId="search-input"
      value={value}
      onValueChange={onChange}
      {...(ariaLabel !== undefined ? { "aria-label": ariaLabel } : {})}
      {...(placeholder !== undefined ? { placeholder } : {})}
      className={className ?? ""}
      leading={<MagnifierSvg />}
      inputRef={inputRef}
      {...(value !== "" ? { trailing: clearButton } : {})}
    />
  );
}
