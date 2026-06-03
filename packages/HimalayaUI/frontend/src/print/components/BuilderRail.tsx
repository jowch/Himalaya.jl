import type { ReactNode } from "react";
import { Button, IconButton, Kicker, SegmentedControl, Slider, Field } from "../ui";
import { AutoGroup } from "./AutoGroup";
import { RailSection } from "./RailSection";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

type Scale = "log" | "lin";

export interface BuilderRailProps {
  grouping: ReactNode;
  onConfirm?: () => void;
  onAdjust?: () => void;
  orderedBy: string;
  orderNote?: ReactNode;
  onChangeOrder?: () => void;
  offset: number;
  onOffsetChange: (v: number) => void;
  scale: Scale;
  onScaleChange: (s: Scale) => void;
  traces: ReactNode;
  onAddSample?: () => void;
  onCopyPng?: () => void;
  onCollapse?: () => void;
  /** PLACEMENT-ONLY. Width is page-owned (the assembly sets it). */
  className?: string;
}

/**
 * BuilderRail — the series-builder "Compose" editing rail (mockup `.rail`).
 *
 * Presentational contract (Batch 7/9/10/11/12): holds NO local state. offset /
 * scale / the trace rows (the `traces` slot) / every handler are PROPS; the
 * SeriesBuilderAssembly page owns the real state and the figure↔rail link.
 *
 * LOAD-BEARING OMISSIONS ("controls don't lie"): the mockup's Representation
 * section (Waterfall/Heatmap) and the "Track reflections" toggle drive renderers
 * that are out-of-scope / deferred (HeatmapChart out-of-scope; TrackingLine
 * deferred — same call as Batch 9 SeriesPlate). Both are OMITTED; the Display
 * section keeps only the offset slider + the log/linear scale toggle.
 */
export function BuilderRail({
  grouping,
  onConfirm,
  onAdjust,
  orderedBy,
  orderNote,
  onChangeOrder,
  offset,
  onOffsetChange,
  scale,
  onScaleChange,
  traces,
  onAddSample,
  onCopyPng,
  onCollapse,
  className,
}: BuilderRailProps): JSX.Element {
  return (
    <aside
      data-testid="builder-rail"
      className={cx(
        "flex flex-col gap-5 bg-paper-sunk border-l border-hair px-5 pt-4 pb-7 overflow-y-auto",
        className,
      )}
    >
      <div className="flex items-center justify-between">
        <Kicker tone="faint">Compose</Kicker>
        <IconButton label="Collapse rail" tone="ghost" onClick={onCollapse}>
          &#8250;
        </IconButton>
      </div>

      <AutoGroup
        variant="compose"
        title="Auto-grouped"
        actions={[
          { label: "Confirm series", ...(onConfirm ? { onClick: onConfirm } : {}) },
          { label: "Adjust", muted: true, ...(onAdjust ? { onClick: onAdjust } : {}) },
        ]}
      >
        {grouping}
      </AutoGroup>

      <RailSection label="Ordering variable" {...(orderNote != null ? { note: orderNote } : {})}>
        <Field value={orderedBy} {...(onChangeOrder ? { onClick: onChangeOrder } : {})} />
      </RailSection>

      <RailSection label="Display">
        <Slider
          label="Trace offset"
          valueDisplay={`${offset.toFixed(2)}×`}
          value={offset}
          min={0.4}
          max={1.4}
          step={0.05}
          onChange={onOffsetChange}
        />
        <SegmentedControl
          aria-label="q scale"
          stretch
          options={[
            { value: "log", label: "log q" },
            { value: "lin", label: "linear q" },
          ]}
          value={scale}
          onChange={onScaleChange}
        />
      </RailSection>

      <RailSection label="Traces — drag to reorder">
        <div className="flex flex-col gap-0.5">{traces}</div>
      </RailSection>

      <div className="flex gap-2">
        <Button variant="outline" className="flex-1" {...(onAddSample ? { onClick: onAddSample } : {})}>
          + Add sample
        </Button>
        <Button variant="outline" className="flex-1" {...(onCopyPng ? { onClick: onCopyPng } : {})}>
          Copy as PNG
        </Button>
      </div>

      <div className="text-caption text-ink-faint leading-relaxed">
        The plate above is the figure as it will export. What you compose is what you publish.
      </div>
    </aside>
  );
}
