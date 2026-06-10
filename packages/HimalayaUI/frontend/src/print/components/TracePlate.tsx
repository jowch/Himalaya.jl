import { useEffect, useState, type ReactNode } from "react";
import { Card, Button, HintText, Input, SegmentedControl } from "../ui";
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
  /** Incoming hot q from another panel (q-link sink) → forwarded to TracePlot. */
  hoveredQ?: number;
  /** Emitted when the user hovers a peak (q-link source) → forwarded to TracePlot. */
  onHoverQ?: (q: number | undefined) => void;
  /** When non-empty, peaks NOT in this set dim to neutral → forwarded to TracePlot. */
  highlightPeakIds?: ReadonlySet<number>;
  /** Plot height in px. Default 360. */
  plotHeight?: number;
  /** Extra toolbar actions rendered after the built-in controls (e.g. figure export). */
  actions?: ReactNode;
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
  hoveredQ,
  onHoverQ,
  highlightPeakIds,
  plotHeight = 360,
  actions,
  className,
}: TracePlateProps): JSX.Element {
  // "+ Peak" arm gate: TracePlot edits peaks on a plot click — empty space
  // adds (onAddPeak), a peak removes / alt-excludes (onClickPeak). BOTH are
  // destructive-ish, so the arm governs all peak editing: while disarmed the
  // plot is read-only (hover/q-link still work) and a stray click neither adds
  // NOR removes a peak. Rebuilt explicitly (not `{...interaction, onAddPeak:
  // undefined}`) to satisfy exactOptionalPropertyTypes.
  const gatedInteraction: TracePlotInteraction | false | undefined = !interaction
    ? interaction
    : {
        onXDomain: interaction.onXDomain,
        ...(addPeakArmed && interaction.onAddPeak
          ? { onAddPeak: interaction.onAddPeak }
          : {}),
        ...(addPeakArmed && interaction.onClickPeak
          ? { onClickPeak: interaction.onClickPeak }
          : {}),
        ...(interaction.onReset ? { onReset: interaction.onReset } : {}),
        ...(interaction.hitTolerancePx !== undefined
          ? { hitTolerancePx: interaction.hitTolerancePx }
          : {}),
      };

  // ── armed add-at-q (keyboard parity for click-empty-space-adds) ────────────
  // Validate against the trace's q extent: a peak outside the measured window
  // would be unanchorable, so Add stays disabled there.
  const [qText, setQText] = useState("");
  let qMin = Infinity;
  let qMax = -Infinity;
  for (const q of trace.trace.q) {
    if (!Number.isFinite(q)) continue;
    if (q < qMin) qMin = q;
    if (q > qMax) qMax = q;
  }
  const qParsed = Number(qText);
  const qValid =
    qText.trim() !== "" &&
    Number.isFinite(qParsed) &&
    qParsed >= qMin &&
    qParsed <= qMax;
  const onAddPeakAtQ =
    addPeakArmed && interaction ? interaction.onAddPeak : undefined;

  // A typed-but-unsubmitted q must not survive a disarm/re-arm round trip —
  // it would reappear as stale state. Clear it whenever the plate disarms.
  useEffect(() => {
    if (!addPeakArmed) setQText("");
  }, [addPeakArmed]);

  // The hint promises only the verbs that are actually wired: with no
  // onAddPeak the add sentence would be false, and with no onClickPeak the
  // remove sentence would be worse than false (TracePlot's click handler
  // falls through to onAddPeak on a peak hit, so "remove" would duplicate-add).
  const hintSentences = interaction
    ? [
        ...(interaction.onAddPeak ? ["Click the trace to add a peak."] : []),
        ...(interaction.onClickPeak
          ? ["Click a peak to remove it.", "Alt-click excludes it from indexing."]
          : []),
      ]
    : [];

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
          {actions}
        </ToolBar>
      </PlateHeader>
      {/* Armed-only edit legend: only the WIRED pointer verbs, plus the
          keyboard add path. Hidden while disarmed — the plot is read-only
          then and a standing instruction would lie. */}
      {addPeakArmed && hintSentences.length > 0 && (
        <div
          data-testid="peak-edit-hint"
          className="mt-2 flex flex-wrap items-center justify-between gap-x-4 gap-y-1.5"
        >
          <HintText className="min-w-0">{hintSentences.join(" ")}</HintText>
          {onAddPeakAtQ && (
            <form
              data-testid="add-peak-at-q"
              className="flex shrink-0 items-center gap-1.5"
              onSubmit={(e) => {
                e.preventDefault();
                if (!qValid) return;
                onAddPeakAtQ(qParsed);
                setQText("");
              }}
            >
              <Input
                value={qText}
                onValueChange={setQText}
                type="number"
                inputMode="decimal"
                step="any"
                min={qMin}
                max={qMax}
                mono
                inputSize="sm"
                invalid={qText.trim() !== "" && !qValid}
                aria-label="q value for new peak"
                placeholder="q"
                testId="add-peak-q-input"
                className="w-24"
              />
              <Button
                type="submit"
                variant="outline"
                disabled={!qValid}
                aria-label="Add peak at q"
              >
                Add
              </Button>
            </form>
          )}
        </div>
      )}
      <TracePlot
        trace={trace}
        height={plotHeight}
        xType={scale === "log" ? "log" : "linear"}
        axes
        {...(xDomain !== undefined ? { xDomain } : {})}
        {...(gatedInteraction !== undefined ? { interaction: gatedInteraction } : {})}
        {...(hoveredQ !== undefined ? { hoveredQ } : {})}
        {...(onHoverQ !== undefined ? { onHoverQ } : {})}
        {...(highlightPeakIds !== undefined ? { highlightPeakIds } : {})}
        paperColor="var(--color-plate)"
        data-testid="trace-plate-plot"
        className="mt-2"
      />
    </Card>
  );
}
