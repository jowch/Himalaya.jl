import { Kicker } from "./ui/Kicker";

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
      <Kicker tone="accent" data-testid="focus-plot-kicker">Integration</Kicker>
      <h1
        data-testid="focus-plot-title"
        className="text-display leading-tight text-ink truncate max-w-[44ch]"
      >
        {sampleName}
      </h1>
      {/* R3-N2 (#209): subline previously capped at `max-w-[60ch]` + `.truncate`,
          which cut the most diagnostic suffix ("representative exposure …")
          off mid-token on real samples. Raise to 80ch — the suffix carries the
          actual filename stem (the load-bearing provenance bit), so a truncate
          at 60ch hides the highest-signal segment. 80ch comfortably fits
          smp_NN + beamtime + "representative exposure smp_NN_eNN" without
          wrapping below `xl`. */}
      <div
        data-testid="focus-plot-sub"
        className="mt-1 font-mono text-xs text-ink-faint truncate max-w-[80ch]"
      >
        {segments.join(" · ")}
      </div>
    </div>
  );
}
