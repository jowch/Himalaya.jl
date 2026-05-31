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
};

export type CardProps<T extends CardElement = "div"> = CardOwnProps<T> &
  Omit<ElementProps[T], keyof CardOwnProps<T>>;

const borderClass: Record<NonNullable<CardOwnProps<CardElement>["border"]>, string> = {
  hair: "border border-hair",
  strong: "border border-hair-strong",
};

export function Card<T extends CardElement = "div">({
  as,
  elevated = false,
  border = "hair",
  draft = false,
  className = "",
  children,
  ...rest
}: CardProps<T>): JSX.Element {
  const Tag = (as ?? "div") as CardElement;
  // elevated → `.card` (single source of the shadow). flat → tonal hairline,
  // never a shadow. Radius token (rounded-md = 5px) comes from Phase 0.
  const appearance = draft
    ? "rounded-md bg-plate border border-dashed border-hair-strong"
    : elevated
      ? "card"
      : `rounded-md bg-plate ${borderClass[border]}`;
  const props = {
    className: `${appearance} ${className}`.trim(),
    ...(draft ? { "data-draft": "true" } : {}),
    ...(elevated && !draft ? { "data-elevated": "true" } : {}),
    // Default an explicit button type so a Card-as-button inside a <form> never
    // silently acts as a submit; a consumer-passed `type` (via ...rest) still wins.
    ...(Tag === "button" ? { type: "button" } : {}),
    ...rest,
  } as Record<string, unknown>;
  return <Tag {...(props as HTMLAttributes<HTMLElement>)}>{children}</Tag>;
}
