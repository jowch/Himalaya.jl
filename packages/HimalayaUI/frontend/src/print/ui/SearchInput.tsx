interface SearchInputProps {
  value: string;
  onChange: (v: string) => void;
  placeholder?: string;
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** Standard product search field: leading magnifier glyph + transparent text
 *  input on a plate. The accent ring lands on the whole field via
 *  `focus-within:border-accent` (the rationed accent-as-focus-mark, applied to
 *  the field rather than the bare input). The fixed `min-w-[230px]` is an
 *  allowed field dimension in print/ui.
 *
 *  C — focus is the `focus-within:border-accent` ring; hover N/A (a field, not a
 *  button). D — `py-1.5` + min-width row gives a comfortable target. E — input
 *  text is sans/base (prose entry), correct. F — magnifier stroke sourced from a
 *  token (`var(--color-ink-faint)`), never hex. */
export function SearchInput({
  value,
  onChange,
  placeholder,
  className,
}: SearchInputProps): JSX.Element {
  return (
    <div
      data-testid="search-input"
      className={cx(
        "inline-flex items-center gap-2 bg-plate border border-hair-strong rounded-sm px-3 py-1.5 min-w-[230px] transition-colors",
        "focus-within:border-accent",
        className,
      )}
    >
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
      <input
        type="search"
        value={value}
        onChange={(e) => onChange(e.target.value)}
        placeholder={placeholder}
        className="flex-1 bg-transparent border-none outline-none text-base text-ink placeholder:text-ink-faint"
      />
    </div>
  );
}
