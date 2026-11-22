import type {
  ButtonHTMLAttributes,
  HTMLAttributes,
  LiHTMLAttributes,
} from "react";

export type CardElement = "div" | "button" | "section" | "li" | "article";

type ElementProps = {
  div: HTMLAttributes<HTMLDivElement>;
  section: HTMLAttributes<HTMLElement>;
  li: LiHTMLAttributes<HTMLLIElement>;
  button: ButtonHTMLAttributes<HTMLButtonElement>;
  article: HTMLAttributes<HTMLElement>;
};

type CardOwnProps<T extends CardElement> = {
  as?: T;
  elevated?: boolean;
  border?: "hair" | "strong";
  draft?: boolean;
  padding?: "sm" | "md" | "lg";
  selected?: boolean;
  previewed?: boolean;
  interactive?: boolean;
};
export type CardProps<T extends CardElement = "div"> = CardOwnProps<T> &
  Omit<ElementProps[T], keyof CardOwnProps<T> | "title">;
const borderClass: Record<NonNullable<CardOwnProps<CardElement>["border"]>, string> = {
  hair: "border border-hair",
  strong: "border border-hair-strong",
};
const paddingClass: Record<NonNullable<CardOwnProps<CardElement>["padding"]>, string> = {
  sm: "p-3",
  md: "p-4",
  lg: "p-7",
};
export function Card<T extends CardElement = "div">({
  as,
  elevated = false,
  border = "hair",
  draft = false,
  padding,
  selected = false,
  previewed = false,
  interactive = false,
  className = "",
  children,
  ...rest
}: CardProps<T>): JSX.Element {
  const Tag = (as ?? "div") as CardElement;
  const appearance =
    (draft
      ? "rounded-md bg-plate border border-dashed border-hair-strong"
      : elevated
        ? "card"
        : `rounded-md bg-plate ${borderClass[border]}`) +
    (interactive ? " cursor-pointer hover:border-hair-strong" : "") +
    (padding ? ` ${paddingClass[padding]}` : "");
  const ringStyle: { borderColor?: string; boxShadow?: string } = {};
  if (selected) {
    ringStyle.borderColor = "color-mix(in oklab, var(--color-accent) 42%, var(--color-hair))";
    ringStyle.boxShadow = "inset 0 0 0 1px color-mix(in oklab, var(--color-accent) 20%, transparent)";
  }
  if (previewed) {
    const previewRing = "inset 0 0 0 2px color-mix(in oklab, var(--color-ink-soft) 30%, transparent)";
    ringStyle.boxShadow = ringStyle.boxShadow ? `${ringStyle.boxShadow}, ${previewRing}` : previewRing;
  }
  const hasRingStyle = ringStyle.borderColor !== undefined || ringStyle.boxShadow !== undefined;
  const props = {
    className: `${appearance} ${className}`.trim(),
    ...(hasRingStyle ? { style: ringStyle } : {}),
    ...(draft ? { "data-draft": "true" } : {}),
    ...(elevated && !draft ? { "data-elevated": "true" } : {}),
    ...(padding ? { "data-padding": padding } : {}),
    ...(selected ? { "data-selected": "true" } : {}),
    ...(previewed ? { "data-previewed": "true" } : {}),
    ...(interactive ? { "data-interactive": "true" } : {}),
    ...(Tag === "button" ? { type: "button" } : {}),
    ...rest,
  } as Record<string, unknown>;
  return <Tag {...(props as HTMLAttributes<HTMLElement>)}>{children}</Tag>;
}