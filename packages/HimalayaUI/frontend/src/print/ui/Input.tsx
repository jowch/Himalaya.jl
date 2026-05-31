import type { InputHTMLAttributes, ReactNode } from "react";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface InputProps extends Omit<InputHTMLAttributes<HTMLInputElement>, "size"> {
  value: string;
  onValueChange: (v: string) => void;
  inputSize?: "sm" | "md";
  /** Adornment slot, e.g. a search magnifier svg. */
  leading?: ReactNode;
  /** Adornment slot, e.g. a clear button or unit suffix. */
  trailing?: ReactNode;
  invalid?: boolean;
  /** Overrides the wrapper's `data-testid`, default `"input"`; lets SearchInput
   *  etc. keep their contract. */
  testId?: string;
}

const sizeClass: Record<"sm" | "md", string> = {
  sm: "px-2.5 py-1",
  md: "px-3 py-1.5",
};

/** Base recessed text field (DESIGN.md §5 Inputs): a `plate` well with a
 *  `hair-strong` hairline and 5px radius. The accent ring lands on the whole
 *  field via `focus-within:border-accent` (the rationed accent-as-focus-mark,
 *  applied to the field rather than the bare input). The shared base for
 *  SearchInput, the scoping order-field, note entry, and the tag editor.
 *
 *  C — focus is the `focus-within:border-accent` ring; hover N/A (a field, not a
 *  button); disabled flows through `...rest`. Invalid = `border-error`, to be
 *  second-channeled by consumer error text (the field carries `aria-invalid` for
 *  assistive tech). F — appearance from `@theme` tokens (`bg-plate`,
 *  `border-hair-strong`, `border-error`, `text-ink`), no hex, 5px `rounded-sm`. */
export function Input({
  value,
  onValueChange,
  inputSize = "md",
  leading,
  trailing,
  invalid = false,
  testId,
  className = "",
  ...rest
}: InputProps): JSX.Element {
  return (
    <div
      data-testid={testId ?? "input"}
      data-invalid={invalid ? "true" : undefined}
      className={cx(
        "inline-flex items-center gap-2 bg-plate border rounded-sm transition-colors focus-within:border-accent",
        invalid ? "border-error" : "border-hair-strong",
        sizeClass[inputSize],
        className,
      )}
    >
      {leading}
      <input
        value={value}
        onChange={(e) => onValueChange(e.target.value)}
        aria-invalid={invalid || undefined}
        className="flex-1 bg-transparent border-none outline-none text-base text-ink placeholder:text-ink-faint min-w-0"
        {...rest}
      />
      {trailing}
    </div>
  );
}
