import { forwardRef } from "react";
import type { InputHTMLAttributes, MutableRefObject, ReactNode, Ref } from "react";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** Apply one node to several refs (forwarded `ref` + the `inputRef` prop). */
function setRefs<T>(node: T | null, ...refs: Array<Ref<T> | undefined>): void {
  for (const r of refs) {
    if (typeof r === "function") r(node);
    else if (r) (r as MutableRefObject<T | null>).current = node;
  }
}

export interface InputProps extends Omit<InputHTMLAttributes<HTMLInputElement>, "size"> {
  value: string;
  onValueChange: (v: string) => void;
  inputSize?: "sm" | "md";
  /**
   * "field" (default): the recessed plate well below. "title": a plate/figure
   * title that is editable in place — Newsreader serif display ink with NO box,
   * just a dotted underline that firms to the accent on focus. The DESIGN.md
   * "confident ink, re-openable (not a permanent open field)" idiom, shared with
   * the scoping plate's title. Serif-Means-Title: only use for an actual title.
   */
  variant?: "field" | "title";
  /** Adornment slot, e.g. a search magnifier svg. */
  leading?: ReactNode;
  /** Adornment slot, e.g. a clear button or unit suffix. */
  trailing?: ReactNode;
  invalid?: boolean;
  /** Renders the inner input in monospace — for lattice/numeric fields. */
  mono?: boolean;
  /** Overrides the wrapper's `data-testid`, default `"input"`; lets SearchInput
   *  etc. keep their contract. */
  testId?: string;
  /** Ref to the inner <input> — lets composers (e.g. SearchInput's clear
   *  affordance) return focus to the field after acting on it. */
  inputRef?: Ref<HTMLInputElement>;
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
export const Input = forwardRef<HTMLInputElement, InputProps>(function Input({
  value,
  onValueChange,
  inputSize = "md",
  leading,
  trailing,
  invalid = false,
  mono = false,
  variant = "field",
  testId,
  inputRef,
  className = "",
  ...rest
}, ref): JSX.Element {
  const isTitle = variant === "title";
  return (
    <div
      data-testid={testId ?? "input"}
      data-variant={isTitle ? "title" : undefined}
      data-invalid={invalid ? "true" : undefined}
      className={cx(
        "inline-flex items-center gap-2 transition-colors",
        isTitle
          // Title: no well, no box — a dotted underline that firms to the
          // accent on focus (confident ink, re-openable; Serif-Means-Title).
          ? cx("border-b border-dotted pb-px focus-within:border-accent",
              invalid ? "border-error" : "border-hair-strong")
          : cx("bg-plate border rounded-sm focus-within:border-accent",
              invalid ? "border-error" : "border-hair-strong",
              sizeClass[inputSize]),
        className,
      )}
    >
      {leading}
      <input
        ref={(node) => setRefs(node, ref, inputRef)}
        value={value}
        onChange={(e) => onValueChange(e.target.value)}
        aria-invalid={invalid || undefined}
        className={cx(
          "flex-1 bg-transparent border-none outline-none placeholder:text-ink-soft min-w-0",
          isTitle
            // Title font scales with inputSize: sm → headline (a compact inline
            // rename, e.g. the grouping SampleFold), md (default) → display.
            ? cx(inputSize === "sm" ? "text-headline" : "text-display", "text-ink")
            : "text-base text-ink",
          mono && "font-mono",
        )}
        {...rest}
      />
      {trailing}
    </div>
  );
});
