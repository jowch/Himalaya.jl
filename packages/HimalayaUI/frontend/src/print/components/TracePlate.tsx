import type { ReactNode } from "react";
import { Card, Button, SegmentedControl } from "../ui";
import { PlateHeader } from "./PlateHeader";
import { ToolBar } from "./ToolBar";
import { TracePlot, type TraceModel, type TracePlotInteraction } from "../plot/TracePlot";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export type TraceScale = "log" | "lin";

export interface TracePlateProps {
  /** Accent eyebrow, e.g. "Integration". */
  kicker?: ReactNode;
  /** Serif plate title (the sample name). */
  title: ReactNode;
  /** Mono subtitle line (ids · facility · representative exposure). */
  subtitle?: ReactNode;
  /** The fully-built trace model (peaks + phase resolved by the caller). */
  trace: TraceModel;
  /** Current x-scale mode; the consumer owns it. */
  scale: TraceScale;
  onScaleChange: (next: TraceScale) => void;
  /** Auto-fit action. Omit → button hidden. */
  onAutoFit?: () => void;
  /** Add-peak armed (terracotta) toggle state + handler. */
  addPeakArmed?: boolean;
  onToggleAddPeak?: () => void;
  /** Forwarded TracePlot interaction (scroll-zoom / add / select / reset).
   *  Scroll-to-zoom EMITS through `interaction.onXDomain`; for the zoom to
   *  RENDER, the consumer must round-trip the emitted window back via `xDomain`
   *  (TracePlot is fully controlled on the domain — it holds no internal zoom
   *  state). So a zoomable plate is the controlled pair `xDomain` + the
   *  `interaction.onXDomain` handler, both owned by the page. */
  interaction?: TracePlotInteraction | false;
  /** Controlled visible q-window; `null`/omitted = full data extent. Feed the
   *  domain emitted by `interaction.onXDomain` back here to render scroll-zoom. */
  xDomain?: [number, number] | null;
  /** Plot height in px. Default 360. */
  plotHeight?: number;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function TracePlate({
  kicker,
  title,
  subtitle,
  trace,
  scale,
  onScaleChange,
  onAutoFit,
  addPeakArmed = false,
  onToggleAddPeak,
  interaction,
  xDomain,
  plotHeight = 360,
  className,
}: TracePlateProps): JSX.Element {
  return (
    <Card
      as="section"
      elevated
      data-testid="trace-plate"
      className={cx("px-6 pt-5 pb-4", className)}
    >
      <PlateHeader kicker={kicker} title={title} subtitle={subtitle} as="h1">
        <ToolBar>
          <SegmentedControl
            options={[
              { value: "log", label: "log q" },
              { value: "lin", label: "linear q" },
            ]}
            value={scale}
            onChange={onScaleChange}
            aria-label="q scale"
          />
          {onAutoFit && (
            <Button variant="ghost" onClick={onAutoFit}>
              Auto-fit
            </Button>
          )}
          {onToggleAddPeak && (
            <Button armed={addPeakArmed} onClick={onToggleAddPeak}>
              + Peak
            </Button>
          )}
        </ToolBar>
      </PlateHeader>
      <TracePlot
        trace={trace}
        height={plotHeight}
        xType={scale === "log" ? "log" : "linear"}
        axes
        {...(xDomain !== undefined ? { xDomain } : {})}
        {...(interaction !== undefined ? { interaction } : {})}
        paperColor="var(--color-plate)"
        data-testid="trace-plate-plot"
        className="mt-2"
      />
    </Card>
  );
}
