import type {
  ButtonHTMLAttributes,
  HTMLAttributes,
  LiHTMLAttributes,
} from "react";

export type CardElement = "div" | "button" | "section" | "li";

type ElementProps = {
  div: HTMLAttributes<HTMLDivElement>;
  section: HTMLAttributes<HTMLElement>;
  li: LiHTMLAttributes<HTMLLIElement>;
  button: ButtonHTMLAttributes<HTMLButtonElement>;
};

type CardOwnProps<T extends CardElement> = {
  /** Render element. Default "div". The folio card is a <button>; masonry
   *  tiles are <li>/<div>. Placement only — appearance is fixed. */
  as?: T;
  /** true → the single lifted "Plate Lift" object (applies the `.card` CSS
   *  rule: bg-plate + hairline + radius + inset highlight + soft warm drop
   *  shadow). false (default) → flat hairline card: bg-plate + hairline +
   *  radius, NO shadow (Flat-Except-the-Plate). */
  elevated?: boolean;
  /** hairline weight for the FLAT variant. "hair" (default) for inner cards;
   *  "strong" for outer/standalone surfaces. Ignored when elevated (`.card`
   *  owns the border). */
  border?: "hair" | "strong";
  /** Folio draft/recipe card: dashed hair-strong border, no Plate-Lift shadow,
   *  figure dimmed by the consumer. (mockup .card.is-draft) */
  draft?: boolean;
  /** Optional internal padding. Default none (consumers often pad the body themselves; mockup padding varies 13–30px). sm≈13px md≈16px lg≈28px. */
  padding?: "sm" | "md" | "lg";
  /** When true, the card renders the accent selection chrome: accent-tinted
   *  border + a subtle inset accent ring (mirrors the mockup .cand.in signal).
   *  Flat — NO shadow. Composes with any `as` element and with the flat
   *  (non-elevated) variant. Adds `data-selected="true"` for tests and
   *  targeted CSS. */
  selected?: boolean;
};

export type CardProps<T extends CardElement = "div"> = CardOwnProps<T> &
  Omit<ElementProps[T], keyof CardOwnProps<T>>;

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
  className = "",
  children,
  ...rest
}: CardProps<T>): JSX.Element {
  const Tag = (as ?? "div") as CardElement;
  // elevated → `.card` (single source of the shadow). flat → tonal hairline,
  // never a shadow. Radius token (rounded-md = 5px) comes from Phase 0.
  // selected → accent border + inset accent ring (flat, no shadow).
  const appearance =
    (draft
      ? "rounded-md bg-plate border border-dashed border-hair-strong"
      : elevated
        ? "card"
        : `rounded-md bg-plate ${borderClass[border]}`) +
    (padding ? ` ${paddingClass[padding]}` : "");
  // The selected ring style is injected inline so it can use CSS custom
  // properties without minting a new Tailwind class outside print/ui.
  // color-mix is allowed here (this file IS design-guard-exempt: src/print/ui/**).
  const selectedStyle = selected
    ? {
        borderColor: "color-mix(in oklab, var(--color-accent) 42%, var(--color-hair))",
        boxShadow: "inset 0 0 0 1px color-mix(in oklab, var(--color-accent) 20%, transparent)",
      }
    : undefined;
  const props = {
    className: `${appearance} ${className}`.trim(),
    ...(selectedStyle ? { style: selectedStyle } : {}),
    ...(draft ? { "data-draft": "true" } : {}),
    ...(elevated && !draft ? { "data-elevated": "true" } : {}),
    ...(padding ? { "data-padding": padding } : {}),
    ...(selected ? { "data-selected": "true" } : {}),
    // Default an explicit button type so a Card-as-button inside a <form> never
    // silently acts as a submit; a consumer-passed `type` (via ...rest) still wins.
    ...(Tag === "button" ? { type: "button" } : {}),
    ...rest,
  } as Record<string, unknown>;
  return <Tag {...(props as HTMLAttributes<HTMLElement>)}>{children}</Tag>;
}
