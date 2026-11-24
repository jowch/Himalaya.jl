import type { ReactNode } from "react";
import { cx } from "../lib/cx";

type Pressed = boolean | "mixed";

const LOOK = {
  rest: "border-hair-strong bg-plate text-ink-soft",
  onFill: "border-current/40 bg-current/15 text-current",
  pressed: "border-accent-fill bg-accent-fill font-semibold text-on-accent",
  mixed: "border-accent bg-plate font-semibold text-accent",
} as const;

function lookOf(onFill: boolean, pressed: Pressed): string {
  if (pressed === true) return LOOK.pressed;
  if (pressed === "mixed") return LOOK.mixed;
  return onFill ? LOOK.onFill : LOOK.rest;
}

export function Kbd({ children, onFill = false, pressed = false }: { children: ReactNode; onFill?: boolean; pressed?: Pressed }): JSX.Element {
  return (
    <kbd
      data-pressed={pressed === false ? undefined : String(pressed)}
      className={cx(
        "inline-flex min-w-[1.5em] whitespace-nowrap justify-center rounded-sm border border-b-2 px-1 text-xs leading-[18px]",
        lookOf(onFill, pressed),
      )}
    >
      {children}
    </kbd>
  );
}
