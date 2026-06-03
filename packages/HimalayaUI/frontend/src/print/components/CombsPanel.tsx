import type { ReactNode } from "react";
import { Card, SegmentedControl } from "../ui";
import { CombChart, ResidualChart, type CombSeries } from "../comb";
import { PanelHeader } from "./PanelHeader";
import { CombLegend } from "./CombLegend";

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
  /** Section label. Default "Reflections — comb". */
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
  label = "Reflections — comb",
  className,
}: CombsPanelProps): JSX.Element {
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
      <PanelHeader label={label}>
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
          <CombChart assigned={assigned} leftover={leftover} {...hover} />
        ) : (
          <ResidualChart assigned={assigned} {...hover} />
        )}
      </div>
      <CombLegend className="mt-2.5" />
    </Card>
  );
}
