import { Link, useLocation, useNavigate, useSearchParams } from "react-router-dom";
import { useAppState } from "../state";
import { useCorpusSamples, useExperiments } from "../queries";
import { sampleDisplayName } from "../lib/sample/displayName";
import { TopBar } from "../print/ui/TopBar";
import { Wordmark } from "../print/ui/Wordmark";
import { Kicker } from "../print/ui/Kicker";
import { IconButton } from "../print/ui/IconButton";

interface Stage {
  id: "samples" | "series";
  label: string;
  /** Absolute path the tab links to. */
  to: string;
}

// Samples (#160) and Series (#173) are the two navigable workflow surfaces. The
// focus workspace at `/sample/:id` is NOT a top-level stage — it is reached only
// by opening a sample from the contact sheet, so it has no tab (the right-side
// sample stepper carries the on-focus context instead).
const STAGES: readonly Stage[] = [
  { id: "samples", label: "Samples", to: "/samples" },
  { id: "series", label: "Series", to: "/series" },
];

/**
 * CorpusTopbar — the topbar for the redesigned corpus-scoped shell: a
 * corpus-level wordmark, the three workflow stage-tabs, and the Beamtime facet
 * chip. This matches the corpus-app-shell mockup exactly (spec
 * 2026-05-17-corpus-app-shell-design.md:121) — no utility controls.
 *
 * The Beamtime chip is an experiment picker: it reads and writes the
 * `?beamtime=<experiment_id>` URL query that the /samples contact sheet
 * (#160) filters on. The URL is the only channel — no prop coupling.
 *
 * I5.1 (#182): the legacy AppHeader/UtilityCluster (theme toggle + switch-user
 * avatar) are retired with the dual-nav shell, and are NOT re-homed here — the
 * mockup omits them. R0a (#221) then retired the dark theme entirely, so the
 * `T` theme shortcut is gone too ("The Print" is the single identity). The
 * `username` identity state + persist survive untouched. A multiplayer
 * switch-user control is deferred to a follow-up issue.
 */
export function CorpusTopbar(): JSX.Element {
  const [searchParams, setSearchParams] = useSearchParams();
  const { pathname } = useLocation();
  const navigate = useNavigate();
  const experimentsQuery = useExperiments();
  const beamtime = searchParams.get("beamtime") ?? "";

  // ── Focus-surface affordances (R5 / #228) ─────────────────────────────────
  // The focus workspace lives at `/sample/:id`. The topbar is global, so the
  // per-sample stepper only appears there.
  const onSampleRoute = pathname.startsWith("/sample/");
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const corpusQ = useCorpusSamples();

  // F-13 stepper: order the active sample's experiment-siblings by their corpus
  // order (matches the `,`/`.` shortcut's experiment-scoped semantics). The URL
  // is the focus surface's source of truth, so prev/next navigate the route
  // (one-way URL→store sync stays intact — see useSyncActiveSampleFromRoute).
  const activeSample = activeSampleId !== undefined
    ? corpusQ.data?.find((s) => s.id === activeSampleId)
    : undefined;
  const siblings = activeSample !== undefined
    ? (corpusQ.data ?? []).filter((s) => s.experiment_id === activeSample.experiment_id)
    : [];
  const stepIdx = activeSample !== undefined
    ? siblings.findIndex((s) => s.id === activeSample.id)
    : -1;
  const prevSample = stepIdx > 0 ? siblings[stepIdx - 1] : undefined;
  const nextSample = stepIdx >= 0 && stepIdx < siblings.length - 1
    ? siblings[stepIdx + 1] : undefined;
  const showStepper = onSampleRoute && activeSample !== undefined && stepIdx >= 0;

  function handlePick(event: React.ChangeEvent<HTMLSelectElement>): void {
    const value = event.target.value;
    setSearchParams((prev) => {
      const next = new URLSearchParams(prev);
      if (value === "") next.delete("beamtime");
      else next.set("beamtime", value);
      return next;
    });
  }

  // ── Wordmark — a real router Link home to /samples that carries the contract
  // testid; the `Wordmark` print primitive owns the brand appearance (SANS,
  // 700, wide tracking) and renders the faint " · SAXS" tail. Placement +
  // focus-ring are the only call-site classes.
  const wordmark = (
    <Link
      to="/samples"
      data-testid="corpus-wordmark"
      aria-label="Himalaya SAXS, go to the corpus"
      className="rounded focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent focus-visible:outline-offset-2"
    >
      <Wordmark tail="SAXS">Himalaya</Wordmark>
    </Link>
  );

  // ── Stage tabs — router Links with router-derived active state + the accent
  // leading dot. The print StageTabs primitive is a button tablist (state-driven,
  // not route-driven), so it does not fit; these stay Links whose appearance is
  // expressed with design-system TOKEN utilities (text-ink / text-ink-faint /
  // bg-paper-sunk / bg-print-accent / rounded) — tokens, not appearance literals.
  const stageTabs = (
    <nav data-testid="stage-tabs" aria-label="Workflow stages" className="flex gap-0.5">
      {STAGES.map((s) => {
        // Active = this tab's path is the current route's prefix (router-derived).
        const isActive = pathname.startsWith(s.to);
        return (
          <Link
            key={s.id}
            to={s.to}
            data-testid={`stage-tab-${s.id}`}
            data-active={isActive ? "true" : undefined}
            aria-current={isActive ? "page" : undefined}
            className={
              "px-2.5 py-1.5 rounded text-xs font-semibold uppercase " +
              "tracking-wide no-underline " +
              (isActive ? "text-ink bg-paper-sunk" : "text-ink-faint")
            }
          >
            <span
              aria-hidden="true"
              className="inline-block w-1 h-1 rounded-full mr-1.5 align-middle bg-print-accent"
            />
            {s.label}
          </Link>
        );
      })}
    </nav>
  );

  // ── Beamtime facet chip — honored only on the samples surface (it filters the
  // contact sheet, and the loupe back-link preserves it). On /series and
  // /sample/:id it was changeable but inert — hide it there rather than present
  // a control that does nothing.
  const beamtimeChip = pathname.startsWith("/samples") ? (
    <select
      data-testid="beamtime-chip"
      aria-label="Filter to a beamtime"
      value={beamtime}
      onChange={handlePick}
      className="rounded-full border border-hair-strong bg-plate px-2.5 py-1
                 text-xs font-semibold text-ink"
    >
      <option value="">Beamtime, all experiments</option>
      {(experimentsQuery.data ?? []).map((exp) => (
        <option key={exp.id} value={exp.id}>
          {exp.name ?? `Experiment ${exp.id}`}
        </option>
      ))}
    </select>
  ) : null;

  // ── F-13: per-sample stepper — the focus surface's primary inter-sample nav,
  // shown in the TopBar rightSlot. URL-routed so the one-way route→store sync
  // stays intact.
  const stepper = showStepper ? (
    <div
      data-testid="sample-stepper"
      className="flex items-center gap-2 text-ink"
    >
      <IconButton
        label="Previous sample"
        tone="ghost"
        disabled={prevSample === undefined}
        onClick={() => prevSample && navigate(`/sample/${prevSample.id}`)}
        data-testid="sample-stepper-prev"
      >{"‹"}</IconButton>
      <span className="flex flex-col items-end leading-tight">
        <span className="text-xs font-semibold text-ink">
          {sampleDisplayName(activeSample!)}
        </span>
        <Kicker as="span" tone="faint">sample {stepIdx + 1} of {siblings.length}</Kicker>
      </span>
      <IconButton
        label="Next sample"
        tone="ghost"
        disabled={nextSample === undefined}
        onClick={() => nextSample && navigate(`/sample/${nextSample.id}`)}
        data-testid="sample-stepper-next"
      >{"›"}</IconButton>
    </div>
  ) : null;

  // Compose from the print TopBar slot-shell (wordmark → children → spacer →
  // rightSlot). The shell owns the bar appearance (h-14, hairline, paper, the
  // 24px gutters, the flex-1 spacer); the corpus-topbar contract testid is
  // carried via TopBar's data-testid passthrough.
  return (
    <TopBar
      data-testid="corpus-topbar"
      wordmark={wordmark}
      rightSlot={stepper}
    >
      {stageTabs}
      {beamtimeChip}
    </TopBar>
  );
}
