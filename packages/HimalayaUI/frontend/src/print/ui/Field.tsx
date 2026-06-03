export interface FieldProps {
  value: string;
  onClick?: () => void;
  placeholder?: string;
  className?: string; // PLACEMENT ONLY
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** The `.field` ordering-variable control: a bordered, clickable row showing a
 *  value on the left and a `▾` chevron (ink-faint) on the right. It is a button
 *  whose click opens a menu owned by the page (not built here). When `value` is
 *  empty and a `placeholder` is given, the placeholder reads in ink-faint. */
export function Field({ value, onClick, placeholder, className }: FieldProps): JSX.Element {
  const showPlaceholder = value === "" && placeholder != null;
  return (
    <button
      type="button"
      data-testid="field"
      {...(onClick ? { onClick } : {})}
      className={cx(
        "w-full flex items-center justify-between border border-hair-strong bg-plate rounded px-3 py-2 text-meta font-semibold cursor-pointer",
        className,
      )}
    >
      <span className={showPlaceholder ? "text-ink-faint" : "text-ink"}>
        {showPlaceholder ? placeholder : value}
      </span>
      <span className="text-ink-faint" aria-hidden="true">
        ▾
      </span>
    </button>
  );
}
