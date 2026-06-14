import { useCallback, useEffect, useRef, useState, type ReactNode } from "react";
import { Card, Button, HintText, Input, SegmentedControl } from "../ui";
import { PlateHeader } from "./PlateHeader";
import { ToolBar } from "./ToolBar";
import { TracePlot, type TraceModel, type TracePlotInteraction } from "../plot/TracePlot";
import type { PeakFocusRequest } from "../plot/marks/PlotPeaks";

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
  /** One-shot keyboard-focus re-anchor after a destructive peak edit (WCAG
   *  2.4.3) → forwarded to TracePlot. The consumer computes the surviving peak
   *  id at activation time and bumps the nonce. */
  focusRequest?: PeakFocusRequest;
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
  focusRequest,
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
  // Validate against the q extent the user can act on: a peak outside the
  // measured trace would be unanchorable, and one outside the visible zoom
  // window would be invisible (FO-ZOOMEDIT) — Add stays disabled in both cases.
  // Focus-re-anchor fallback target: when a destructive peak edit leaves no
  // surviving mark to take focus, the plot layer calls onFocusFallback and we
  // park focus on the "+ Peak" button (the keyboard user's last stable handle).
  const addPeakButtonRef = useRef<HTMLButtonElement | null>(null);
  const onFocusFallback = useCallback(() => {
    addPeakButtonRef.current?.focus();
  }, []);

  const [qText, setQText] = useState("");
  let qMin = Infinity;
  let qMax = -Infinity;
  for (const q of trace.trace.q) {
    if (!Number.isFinite(q)) continue;
    if (q < qMin) qMin = q;
    if (q > qMax) qMax = q;
  }
  // FO-ZOOMEDIT: validate the add against the VISIBLE window, not the full
  // trace extent. When zoomed (`xDomain` is a controlled sub-range), a q that
  // is in the data but outside the window would add a peak you cannot see until
  // you zoom out. Clamp the window into the data extent so the bound is never
  // wider than the measured trace; with no zoom (`xDomain` null/omitted) it is
  // exactly [qMin, qMax] — unchanged.
  const addLo = xDomain ? Math.max(qMin, Math.min(xDomain[0], xDomain[1])) : qMin;
  const addHi = xDomain ? Math.min(qMax, Math.max(xDomain[0], xDomain[1])) : qMax;
  const qParsed = Number(qText);
  const qValid =
    qText.trim() !== "" &&
    Number.isFinite(qParsed) &&
    qParsed >= addLo &&
    qParsed <= addHi;
  const onAddPeakAtQ =
    addPeakArmed && interaction ? interaction.onAddPeak : undefined;

  // A typed-but-unsubmitted q must not survive a disarm/re-arm round trip —
  // it would reappear as stale state. Clear it whenever the plate disarms.
  useEffect(() => {
    if (!addPeakArmed) setQText("");
  }, [addPeakArmed]);

  // ── F7: Escape disarms the armed mode (the keyboard exit) ──────────────────
  // Precedence: an OPEN MODAL DIALOG wins — ModalShell owns Escape-to-close
  // and stamps preventDefault on the press it consumes, so the closing press
  // arrives here already-defaultPrevented; only a later Escape disarms. The
  // defaultPrevented check is the load-bearing one: in a real browser a
  // microtask checkpoint runs between the document-level close listener and
  // this window-level one, so React has already unmounted the dialog by the
  // time this fires — DOM presence alone misses the closing press (jsdom's
  // synchronous dispatch hides this). The DOM check stays for open dialogs
  // that keep Escape inert (closeOnEsc=false, parent-owned Escape). The
  // add-at-q Input has no local Escape behavior, so Escape there falls
  // through to the disarm (whose effect above also clears the draft q).
  // suppressGlobalKeys is deliberately NOT used here: it suppresses ALL
  // typing contexts, which would make Escape inert exactly where the armed
  // mode parks the focus.
  useEffect(() => {
    if (!addPeakArmed || !onToggleAddPeak) return;
    const onKey = (e: KeyboardEvent): void => {
      if (e.key !== "Escape") return;
      if (e.defaultPrevented) return;
      if (document.querySelector('[role="dialog"][aria-modal="true"]') !== null) return;
      // WCAG 2.4.3 (FO-FOCUSRETURN): disarming strips every peak mark's
      // tabIndex/role and unmounts the add-at-q field, so an Escape exit while
      // one of those armed-only controls holds focus would drop focus to <body>.
      // Re-anchor to the "+ Peak" button — the keyboard user's stable handle —
      // in that case. Focus already on the toolbar button, or off the plate, is
      // left alone (no yank).
      const ae = document.activeElement;
      const onVanishingControl =
        ae?.closest('[data-role="plot-peaks"]') != null ||
        ae?.getAttribute("aria-label") === "q value for new peak";
      onToggleAddPeak();
      if (onVanishingControl) addPeakButtonRef.current?.focus();
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [addPeakArmed, onToggleAddPeak]);

  // The hint promises only the verbs that are actually wired: with no
  // onAddPeak the add sentence would be false, and with no onClickPeak the
  // remove sentence would be worse than false (TracePlot's click handler
  // falls through to onAddPeak on a peak hit, so "remove" would duplicate-add).
  // The click verb is provenance-split downstream (a manual peak removes, an
  // auto peak toggles off) — TracePlate cannot see a peak's source, so the
  // sentence names both outcomes generically rather than lying about one.
  // "Esc exits." needs the toggle handler — that is what Escape disarms through.
  const hintSentences = [
    ...(interaction && interaction.onAddPeak
      ? ["Click the trace to add a peak."]
      : []),
    ...(interaction && interaction.onClickPeak
      ? ["Click a peak to remove it; auto peaks toggle off instead."]
      : []),
    ...(onToggleAddPeak ? ["Esc exits."] : []),
  ];

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
          {onToggleAddPeak && (
            <Button ref={addPeakButtonRef} armed={addPeakArmed} onClick={onToggleAddPeak}>
              + Peak
            </Button>
          )}
          {actions}
        </ToolBar>
      </PlateHeader>
      {/* Disarmed discoverability cue (H10): when editing is available but not
          armed, a first-timer has no signal that peaks are interactive. A quiet
          one-liner points at the affordance. Shown ONLY while disarmed and only
          when the edit toggle exists; never alongside the armed legend. Rides
          ink-soft (AA) — it is informational small text, not decorative. */}
      {!addPeakArmed && onToggleAddPeak && (
        <p data-testid="peak-edit-cue" className="mt-2 text-meta text-ink-soft">
          Arm + Peak to edit peaks.
        </p>
      )}
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
                min={addLo}
                max={addHi}
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
        {...(focusRequest !== undefined ? { focusRequest, onFocusFallback } : {})}
        figureLabel="Integration trace: intensity vs q"
        paperColor="var(--color-plate)"
        data-testid="trace-plate-plot"
        className="mt-2"
      />
    </Card>
  );
}
