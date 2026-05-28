/**
 * FocusPlotHeader — the focus-workspace trace-plate header (mockup
 * focus-workspace.html `.plate-head`, DESIGN.md "The Print" R3 / #226).
 *
 * Replaces the legacy experiment-picker title (PlotCard's prop-less default)
 * in the focus context, where the route seeds only the sample. The header is
 * passed to <PlotCard headerSlot={...} /> so the picker affordance — and its
 * `onTitleClick → openNavModal` branch — never renders here.
 *
 * Layout (matches `.plate-head`):
 *   Integration                         ← terracotta uppercase kicker (L-8)
 *   <Sample name>                        ← Newsreader serif 27px/500 (L-7)
 *   smp · beamtime · representative exposure <e>   ← mono subline
 */
export interface FocusPlotHeaderProps {
  /** Serif display title — the sample's human name. */
  sampleName: string;
  /** Mono subline lead — the sample's internal code (e.g. "smp_09"). */
  sampleCode: string | null;
  /** Mono subline — the owning experiment / beamtime name. */
  beamtime: string | null;
  /** Mono subline — the representative (selected) exposure label. */
  exposureLabel: string | null;
}

export function FocusPlotHeader({
  sampleName,
  sampleCode,
  beamtime,
  exposureLabel,
}: FocusPlotHeaderProps): JSX.Element {
  // Build the subline from only the segments we actually have, so a cold
  // corpus (no beamtime / no exposure yet) degrades to just the code or
  // nothing — never a stray separator or dangling "representative exposure".
  const segments: string[] = [];
  if (sampleCode) segments.push(sampleCode);
  if (beamtime) segments.push(beamtime);
  if (exposureLabel) segments.push(`representative exposure ${exposureLabel}`);

  return (
    <div data-testid="focus-plot-header" className="min-w-0 flex flex-col items-start">
      <div
        data-testid="focus-plot-kicker"
        className="text-xs font-semibold uppercase tracking-wide text-print-accent"
      >
        Integration
      </div>
      <h1
        data-testid="focus-plot-title"
        className="text-display leading-tight text-ink truncate max-w-[44ch]"
      >
        {sampleName}
      </h1>
      <div
        data-testid="focus-plot-sub"
        className="mt-1 font-mono text-xs text-ink-faint truncate max-w-[60ch]"
      >
        {segments.join(" · ")}
      </div>
    </div>
  );
}
