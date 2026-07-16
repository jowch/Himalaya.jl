import { useCallback, useRef, type ReactNode } from "react";
import { Card, Button, SegmentedControl, Tooltip } from "../ui";
import { PlateHeader } from "./PlateHeader";
import { ToolBar } from "./ToolBar";
import { TracePlot, type TraceModel, type TracePlotInteraction } from "../plot/TracePlot";
import type { PeakFocusRequest } from "../plot/marks/PlotPeaks";
import { cx } from "../../lib/cx";


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
  /** One-shot keyboard-focus re-anchor after a destructive peak edit (WCAG
   *  2.4.3) → forwarded to TracePlot. The consumer computes the surviving peak
   *  id at activation time and bumps the nonce. */
  focusRequest?: PeakFocusRequest;
  /** Plot height in px. Default 360. */
  plotHeight?: number;
  /** While true, the plot area renders a quiet skeleton block (the header stays
   *  — the sample name/subtitle are meaningful the moment the exposure resolves,
   *  before the trace/peaks fetch lands). Same height as the live plot, so the
   *  load→loaded swap never shifts layout. */
  loading?: boolean;
  /** Forwarded y-headroom: expands the y-domain top so the tallest peak's
   *  marker clears the ceiling instead of clamping onto the trace. */
  yHeadroom?: number;
  /** Render the trace LINE in the neutral (gray) colour regardless of phase;
   *  the peak markers keep their phase colour. Keeps an arbitrary
   *  coexistence-phase hue off the curve. */
  neutralLine?: boolean;
  /** Show the peak-label layer (the per-peak ratio labels, e.g. on candidate
   *  hover). Labels are carried on each PlotPeak's `label`. */
  showPeakLabels?: boolean;
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
  focusRequest,
  plotHeight = 360,
  yHeadroom,
  neutralLine = false,
  showPeakLabels = false,
  actions,
  className,
  loading = false,
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

  // Focus-re-anchor fallback target: when a destructive peak edit leaves no
  // surviving mark to take focus, the plot layer calls onFocusFallback and we
  // park focus on the "+ Peak" button (the keyboard user's last stable handle).
  const addPeakButtonRef = useRef<HTMLButtonElement | null>(null);
  const onFocusFallback = useCallback(() => {
    addPeakButtonRef.current?.focus();
  }, []);

  // The hint promises only the verbs that are actually wired: with no
  // onAddPeak the add sentence would be false, and with no onClickPeak the
  // remove sentence would be worse than false (TracePlot's click handler
  // falls through to onAddPeak on a peak hit, so "remove" would duplicate-add).
  // The click verb is provenance-split downstream (a manual peak removes, an
  // auto peak toggles off) — TracePlate cannot see a peak's source, so the
  // sentence names both outcomes generically rather than lying about one.
  // "Esc exits." is gated on onToggleAddPeak being wired, but Escape itself is
  // handled by the page's escapeLadder — the page owns Escape-disarm.
  const hintSentences = [
    ...(interaction && interaction.onAddPeak
      ? ["Click the trace to add a peak."]
      : []),
    ...(interaction && interaction.onClickPeak
      ? ["Click a peak to remove it; auto peaks toggle off instead."]
      : []),
    ...(onToggleAddPeak ? ["Esc exits."] : []),
  ];
  // The guidance lives in a hover/focus tooltip on the "+ Peak" button instead
  // of a standing block of text that pushed the plot down (it was the same copy
  // in both states). Disarmed, it leads with the arm step; armed, it is just the
  // live verbs. Empty when no verbs are wired (the tooltip is then suppressed).
  const editTip =
    hintSentences.length === 0
      ? ""
      : (addPeakArmed ? "" : "Arm to edit peaks: ") + hintSentences.join(" ");

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
          {/* FO-TOOLBAR-FLAT: split the view controls (scale, Auto-fit) from the
              edit/export controls (+ Peak, Copy) with a hairline so the flat
              equal-gap run reads as two intents. Hidden when it would orphan an
              edge (no edit/export actions present). */}
          {(onToggleAddPeak || actions) && (
            <span
              aria-hidden="true"
              className="self-stretch w-px bg-hair-strong my-0.5 mx-0.5"
            />
          )}
          {onToggleAddPeak &&
            (editTip ? (
              // Guidance rides a hover/focus tooltip on the button — it no longer
              // pushes the plot down with a standing block of text. (Our own
              // instant dark tooltip, not the slow browser-native title.)
              <Tooltip label={editTip} side="bottom" multiline>
                <Button ref={addPeakButtonRef} armed={addPeakArmed} onClick={onToggleAddPeak}>
                  + Peak
                </Button>
              </Tooltip>
            ) : (
              <Button ref={addPeakButtonRef} armed={addPeakArmed} onClick={onToggleAddPeak}>
                + Peak
              </Button>
            ))}
          {actions}
        </ToolBar>
      </PlateHeader>
      {loading ? (
        <div
          data-testid="trace-plate-skeleton"
          aria-hidden="true"
          className="mt-2 rounded-md bg-paper-sunk animate-pulse"
          style={{ height: plotHeight }}
        />
      ) : (
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
          {...(focusRequest !== undefined ? { focusRequest, onFocusFallback } : {})}
          {...(yHeadroom !== undefined ? { yHeadroom } : {})}
          {...(neutralLine ? { neutralLine: true } : {})}
          {...(showPeakLabels ? { show: { labels: true } } : {})}
          figureLabel="Integration trace: intensity vs q"
          paperColor="var(--color-plate)"
          data-testid="trace-plate-plot"
          className="mt-2"
        />
      )}
    </Card>
  );
}
