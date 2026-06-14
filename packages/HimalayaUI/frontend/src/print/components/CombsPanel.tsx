import type { ReactNode } from "react";
import { Card, SegmentedControl } from "../ui";
import { CombChart, ResidualChart, type CombSeries } from "../comb";
import { PanelHeader } from "./PanelHeader";
import { CombLegend } from "./CombLegend";
import { ResidualLegend } from "./ResidualLegend";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export type CombView = "comb" | "resid";

export interface CombsPanelProps {
  /** Assigned phases → comb rows. */
  assigned: CombSeries[];
  /** Leftover (unindexed) observed-peak q-values; comb view only. */
  leftover?: number[];
  /** Which view is showing; the consumer owns it. */
  view: CombView;
  onViewChange: (next: CombView) => void;
  /** Incoming q-link from the trace/detector. */
  hoveredQ?: number;
  onHoverQ?: (q?: number) => void;
  /** FO-COMB-AXIS: the trace's q-domain, so the comb axis spans the same
   *  q-range as the trace (shared ticks + true reflection positions). */
  qDomain?: [number, number];
  /** Section label override. By default the label derives from `view`:
   *  "Reflections · comb" / "Reflections · indexing space". */
  label?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function CombsPanel({
  assigned,
  leftover = [],
  view,
  onViewChange,
  hoveredQ,
  onHoverQ,
  qDomain,
  label,
  className,
}: CombsPanelProps): JSX.Element {
  // Single source: the same `view` that switches the chart also names the section,
  // so flipping to indexing space never keeps the comb's identity in the header.
  const heading =
    label ?? (view === "comb" ? "Reflections · comb" : "Reflections · indexing space");
  const hover = {
    ...(hoveredQ !== undefined ? { hoveredQ } : {}),
    ...(onHoverQ !== undefined ? { onHoverQ } : {}),
  };
  return (
    <Card
      as="section"
      data-testid="combs-panel"
      data-view={view}
      className={cx("flex flex-col p-4", className)}
    >
      <PanelHeader label={heading}>
        <SegmentedControl
          size="xs"
          options={[
            { value: "comb", label: "comb" },
            { value: "resid", label: "indexing space" },
          ]}
          value={view}
          onChange={onViewChange}
          aria-label="comb view"
        />
      </PanelHeader>
      <div data-testid="combs-body" className="flex-1 min-h-0 overflow-hidden">
        {view === "comb" ? (
          <CombChart assigned={assigned} leftover={leftover} {...(qDomain ? { xDomain: qDomain } : {})} {...hover} />
        ) : (
          <ResidualChart assigned={assigned} {...hover} />
        )}
      </div>
      {/* Per-view legend: each view's glyph vocabulary, never the other's. */}
      {view === "comb" ? (
        <CombLegend className="mt-2.5" />
      ) : (
        <ResidualLegend className="mt-2.5" />
      )}
    </Card>
  );
}
