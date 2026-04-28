import type { ReactNode } from "react";

type Slot = "left" | "center" | "right";

interface Props {
  left:   ReactNode;
  center: ReactNode;
  right:  ReactNode;
  /**
   * Stack order on narrow viewports (<1400px).
   * Default: ["center", "right", "left"] — image/plot first, side context after,
   *   chat last (it's a long history; users want primary content first on mobile).
   */
  mobileOrder?: [Slot, Slot, Slot];
  /**
   * Per-slot wrapper className (e.g. min-heights specific to the card content).
   * The grid itself sets layout/overflow; consumers tune content sizing.
   */
  slotClassName?: Partial<Record<Slot, string>>;
  className?: string;
}

const ORDER_CLASS: Record<Slot, Record<number, string>> = {
  left:   { 1: "order-1", 2: "order-2", 3: "order-3" },
  center: { 1: "order-1", 2: "order-2", 3: "order-3" },
  right:  { 1: "order-1", 2: "order-2", 3: "order-3" },
};

/**
 * WorkspaceGrid — shared three-column layout for IndexPage and InspectPage.
 *
 *   < 1400px  — single column stacked. Order controlled by `mobileOrder`.
 *   ≥ 1400px  — three columns: minmax(320px, 22fr) | 56fr | minmax(320px, 22fr)
 *
 * Side cards never go below 320px even at the 1400px floor; the center column
 * absorbs all residual width. Below 1400px we drop straight to a single column —
 * the in-between two-column tier produced cramped side cards (chat at ~280px).
 *
 * Pure layout: no state, no behavior. The grid template + breakpoints live here
 * so the two pages can never drift.
 */
export function WorkspaceGrid({
  left,
  center,
  right,
  mobileOrder = ["center", "right", "left"],
  slotClassName,
  className,
}: Props): JSX.Element {
  const orderOf = (slot: Slot): string => {
    const idx = mobileOrder.indexOf(slot) + 1; // 1-based
    return ORDER_CLASS[slot][idx] ?? "";
  };

  const baseSlot =
    "card overflow-hidden min-[1400px]:order-none min-[1400px]:min-h-0";

  return (
    <div
      data-testid="workspace-grid"
      className={[
        "min-h-0 grid gap-3 flex-1",
        "grid-cols-1",
        "min-[1400px]:grid-cols-[minmax(320px,22fr)_56fr_minmax(320px,22fr)]",
        "h-auto min-[1400px]:flex-1",
        "min-[1400px]:max-h-[min(700px,calc(100dvh-var(--chrome-h)-1.5rem))]",
        className ?? "",
      ].join(" ")}
    >
      <section
        data-slot="left"
        className={`${baseSlot} ${orderOf("left")} ${slotClassName?.left ?? ""}`}
      >
        {left}
      </section>
      <section
        data-slot="center"
        className={`${baseSlot} ${orderOf("center")} ${slotClassName?.center ?? ""}`}
      >
        {center}
      </section>
      <section
        data-slot="right"
        className={`${baseSlot} ${orderOf("right")} ${slotClassName?.right ?? ""}`}
      >
        {right}
      </section>
    </div>
  );
}
