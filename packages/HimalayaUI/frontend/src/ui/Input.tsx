import { forwardRef } from "react";
import type { InputHTMLAttributes, MutableRefObject, ReactNode, Ref } from "react";
import { cx } from "../lib/cx";

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
  variant?: "field" | "title";
  leading?: ReactNode;
  trailing?: ReactNode;
  invalid?: boolean;
  mono?: boolean;
  testId?: string;
  inputRef?: Ref<HTMLInputElement>;
}
const sizeClass: Record<"sm" | "md", string> = {
  sm: "px-2.5 py-1",
  md: "px-3 py-1.5",
};
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