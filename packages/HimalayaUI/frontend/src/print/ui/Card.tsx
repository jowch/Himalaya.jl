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
  /** Keyboard/hover "cursor" — the candidate the arrows are pointing at (Focus
   *  ↑/↓ preview). A neutral inset ink ring, deliberately DISTINCT from the
   *  accent `selected` ring (previewed ≠ in-the-call). Composes with `selected`
   *  (both rings show). Adds `data-previewed="true"`. */
  previewed?: boolean;
  /** Clickable card / door. The quiet house hover affordance — hairline
   *  firming (`hair` → `hair-strong`) + pointer cursor, riding the global
   *  120ms colour transition. NO motion, NO shadow change. Adds
   *  `data-interactive="true"` for tests. */
  interactive?: boolean;
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
  previewed = false,
  interactive = false,
  className = "",
  children,
  ...rest
}: CardProps<T>): JSX.Element {
  const Tag = (as ?? "div") as CardElement;
  // elevated → `.card` (single source of the Plate-Lift shadow, sourced from
  // the `--shadow-plate-lift` token). flat → tonal hairline, never a shadow.
  // Radius token (rounded-md = 5px) comes from Phase 0.
  // selected → accent border + inset accent ring (flat, no shadow).
  // interactive → quiet hover hairline firming (utilities layer beats the
  // `.card` base-layer border colour), pointer cursor. No motion.
  const appearance =
    (draft
      ? "rounded-md bg-plate border border-dashed border-hair-strong"
      : elevated
        ? "card"
        : `rounded-md bg-plate ${borderClass[border]}`) +
    (interactive ? " cursor-pointer hover:border-hair-strong" : "") +
    (padding ? ` ${paddingClass[padding]}` : "");
  // The selected ring style is injected inline so it can use CSS custom
  // properties without minting a new Tailwind class outside print/ui.
  // color-mix is allowed here (this file IS design-guard-exempt: src/print/ui/**).
  // Selected (accent) and previewed (neutral ink) rings compose — a candidate can
  // be both in-the-call AND under the arrow cursor. color-mix is allowed here
  // (this file IS design-guard-exempt: src/print/ui/**).
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
    // Default an explicit button type so a Card-as-button inside a <form> never
    // silently acts as a submit; a consumer-passed `type` (via ...rest) still wins.
    ...(Tag === "button" ? { type: "button" } : {}),
    ...rest,
  } as Record<string, unknown>;
  return <Tag {...(props as HTMLAttributes<HTMLElement>)}>{children}</Tag>;
}
