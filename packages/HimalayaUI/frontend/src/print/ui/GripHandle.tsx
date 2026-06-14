function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

interface GripHandleProps {
  className?: string;
}

// The drag handle (mockup `.grip`): two ⋮ glyphs packed tight, a quiet STANDING
// reorderable signal at rest that firms on ROW hover. The reveal is driven by
// the consumer row's `group` (we expose it via `group-hover:`); the consumer
// owns the actual drag wiring.
//
// `tracking-tighter` approximates the mockup's −2px column packing on-scale (no
// arbitrary `tracking-[-2px]`).
//
// CC-GRIP-FAINT: the rest tone is `ink-faint` (3.16:1 on paper — clears the 3:1
// non-text-UI bar) so the grip is legible as a standing affordance rather than a
// near-invisible `hair-strong` smudge; it firms to `ink-soft` on hover.
//
// `aria-hidden` is correct here because pointer-drag is not the keyboard-primary
// affordance: the GRIP is purely visual. (Impeccable: G = only a `group-hover`
// COLOR change, no layout/size animation; the color shift is tonal
// ink-faint→ink-soft, not chromatic, so no second-channel concern; C/D — the
// interactive drag/reorder lives in the consumer, N/A at the primitive level.)
//
// CONSUMER RESPONSIBILITY: because this grip is aria-hidden, the consumer row
// MUST provide a keyboard reorder path so reordering is not pointer-only.
// ScopeSampleRow's `onMoveBy` contract is the canonical discharge: when given,
// the row wraps this glyph in a real "Reorder {name}" button handling
// ArrowUp/ArrowDown (SC-KBD). Other reorderable consumers (e.g. the builder's
// MemberRow, tracked as BU-ORDERINERT) should adopt the same contract.
export function GripHandle({ className }: GripHandleProps): JSX.Element {
  return (
    <span
      data-testid="grip-handle"
      aria-hidden="true"
      className={cx(
        "select-none leading-none text-base text-ink-faint group-hover:text-ink-soft cursor-grab flex-shrink-0 tracking-tighter",
        className,
      )}
    >
      ⋮⋮
    </span>
  );
}
