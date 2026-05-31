function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

interface GripHandleProps {
  className?: string;
}

// The drag handle (mockup `.grip`): two ⋮ glyphs packed tight, dim at rest and
// brightening on ROW hover. The reveal is driven by the consumer row's `group`
// (we expose it via `group-hover:`); the consumer owns the actual drag wiring.
//
// `tracking-tighter` approximates the mockup's −2px column packing on-scale (no
// arbitrary `tracking-[-2px]`).
//
// `aria-hidden` is correct here because pointer-drag is not the keyboard-primary
// affordance: the GRIP is purely visual. (Impeccable: G = only a `group-hover`
// COLOR change, no layout/size animation; the color shift is tonal
// hair-strong→ink-faint, not chromatic, so no second-channel concern; C/D — the
// interactive drag/reorder lives in the consumer, N/A at the primitive level.)
//
// PHASE-2 CONSUMER RESPONSIBILITY: because this grip is aria-hidden, the consumer
// row MUST provide a keyboard reorder path (e.g. move-up/down buttons or
// aria-grabbed semantics) so reordering is not pointer-only.
export function GripHandle({ className }: GripHandleProps): JSX.Element {
  return (
    <span
      data-testid="grip-handle"
      aria-hidden="true"
      className={cx(
        "select-none leading-none text-base text-hair-strong group-hover:text-ink-faint cursor-grab flex-shrink-0 tracking-tighter",
        className,
      )}
    >
      ⋮⋮
    </span>
  );
}
