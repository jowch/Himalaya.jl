import { Link, useLocation, useMatch, useNavigate, useSearchParams } from "react-router-dom";
import { useCorpusSamples, useExperiments } from "../../queries";
import { resolveRouteSampleStatus } from "../../hooks/useSyncActiveSampleFromRoute";
import { useExperimentSiblings } from "../../hooks/useExperimentSiblings";
import { sampleDisplayName } from "../../lib/sample/displayName";
import { resolveSampleOrder, sampleNeighbors } from "../../lib/sample/sampleOrder";
import { SampleStepper } from "./SampleStepper";
import { SHORTCUTS } from "./shortcuts";
import {
  resolveExperimentFilter,
  UNKNOWN_BEAMTIME_LABEL,
} from "../../lib/experimentFilter";
import { TopBar } from "../ui/TopBar";
import { Wordmark } from "../ui/Wordmark";

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
  const location = useLocation();
  const { pathname } = location;
  const navigate = useNavigate();
  const experimentsQuery = useExperiments();

  // ── Focus-surface affordances (R5 / #228) ─────────────────────────────────
  // The focus workspace lives at `/sample/:id`. The topbar is global, so the
  // per-sample stepper only appears there.
  const onSampleRoute = pathname.startsWith("/sample/");
  const corpusQ = useCorpusSamples();

  // F-STALEURL honesty gate: the topbar renders ABOVE the routed element, so
  // useParams cannot see the child route's :sampleId; useMatch reads it off
  // the location instead. A bogus param never seeds the store, so the stale
  // activeSampleId would otherwise keep the previous sample's stepper alive
  // over a "Sample not found" body. Hide on "unknown" only: "pending" cannot
  // be judged yet, and gating on param equality with the store would flicker
  // the stepper on every valid step (the store lags the URL by one effect
  // tick, while a known param stays "found" throughout).
  const sampleMatch = useMatch("/sample/:sampleId");
  const routeStatus = resolveRouteSampleStatus(
    sampleMatch?.params.sampleId,
    corpusQ.data,
  );

  // F-13 stepper: the active sample's experiment-siblings in corpus order,
  // via the SHARED useExperimentSiblings derivation (F5) — the `,`/`.` global
  // shortcut steps through the identical list, so the two can never disagree.
  // The URL is the focus surface's source of truth, so prev/next navigate the
  // route (one-way URL→store sync stays intact — see useSyncActiveSampleFromRoute).
  const {
    activeSample,
    siblings,
    index: stepIdx,
    prev: prevSample,
    next: nextSample,
  } = useExperimentSiblings();
  const showStepper =
    onSampleRoute &&
    routeStatus !== "unknown" &&
    activeSample !== undefined &&
    stepIdx >= 0;

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
  // expressed with design-system TOKEN utilities (text-ink / text-ink-soft /
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
              (isActive ? "text-ink bg-paper-sunk" : "text-ink-soft")
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
  //
  // SA-F5 honesty: the SAME resolver the contact sheet uses judges the URL's
  // filter, so the select and the page body can never disagree. A ?beamtime
  // that matches no <option> would otherwise make the controlled select
  // display "all experiments" while the sheet is filtered — surface the actual
  // state as a disabled option instead (the only way OUT of it is picking a
  // real entry, mirroring the sheet's clear-filter CTA).
  const filter = resolveExperimentFilter(
    searchParams.get("beamtime"),
    experimentsQuery.data,
    corpusQ.data,
  );
  const beamtimeUnlisted =
    filter.id !== undefined &&
    !(experimentsQuery.data ?? []).some((e) => e.id === filter.id);
  const beamtimeChip = pathname.startsWith("/samples") ? (
    <select
      data-testid="beamtime-chip"
      aria-label="Filter to a beamtime"
      // The resolver's normalized id, not the raw param: a non-canonical
      // numeral ("007") filters the sheet but matches no option value, so the
      // raw string would display "all experiments" over a filtered corpus.
      value={filter.id !== undefined ? String(filter.id) : ""}
      onChange={handlePick}
      // TOP-CHIPMOBILE: below the tablet breakpoint the wordmark + stage tabs
      // push this chip's start to ~x=253, so its natural ~177px width spills past
      // a 390px viewport. Cap + truncate ONLY below `sm` (min-w-0 lets it shrink
      // in the flex row); at sm+ it is uncapped so full experiment names show on
      // the desktop target. Mobile is not the target viewport (desktop tool).
      className="min-w-0 max-w-28 truncate sm:max-w-none
                 rounded-full border border-hair-strong bg-plate px-2.5 py-1
                 text-xs font-semibold text-ink"
    >
      <option value="">Beamtime, all experiments</option>
      {beamtimeUnlisted && (
        <option value={String(filter.id)} disabled>
          {filter.unknown
            ? UNKNOWN_BEAMTIME_LABEL
            : (filter.name ?? `Experiment ${filter.id}`)}
        </option>
      )}
      {(experimentsQuery.data ?? []).map((exp) => (
        <option key={exp.id} value={exp.id}>
          {exp.name ?? `Experiment ${exp.id}`}
        </option>
      ))}
    </select>
  ) : null;

  const prevKey = SHORTCUTS.prevSample.keys[0]!;
  const nextKey = SHORTCUTS.nextSample.keys[0]!;

  // ── F-13: per-sample stepper — the focus surface's primary inter-sample nav,
  // shown in the TopBar rightSlot. URL-routed so the one-way route→store sync
  // stays intact. Steps the active sample's experiment-siblings.
  const focusStepper = showStepper ? (
    <SampleStepper
      name={sampleDisplayName(activeSample!)}
      index={stepIdx}
      total={siblings.length}
      {...(prevSample ? { onPrev: () => navigate(`/sample/${prevSample.id}`) } : {})}
      {...(nextSample ? { onNext: () => navigate(`/sample/${nextSample.id}`) } : {})}
      prevKey={prevKey}
      nextKey={nextKey}
    />
  ) : null;

  // The SAME stepper on the Loupe (same component + location), stepping the
  // contact-sheet order the loupe was opened with (shared resolveSampleOrder, so
  // it never disagrees with the loupe page's own `[`/`]` keyboard nav).
  const loupeMatch = useMatch("/samples/loupe/:sampleId");
  const loupeSampleId = loupeMatch ? Number(loupeMatch.params.sampleId) : undefined;
  const beamtimeParam = searchParams.get("beamtime");
  const beamtimeNum =
    beamtimeParam !== null && /^\d+$/.test(beamtimeParam) ? Number(beamtimeParam) : undefined;
  const loupeStepper = (() => {
    if (loupeSampleId === undefined || Number.isNaN(loupeSampleId)) return null;
    const corpus = corpusQ.data ?? [];
    const sample = corpus.find((s) => s.id === loupeSampleId);
    if (sample === undefined) return null;
    const sampleOrderState = (location.state as { sampleOrder?: number[] } | null)?.sampleOrder;
    const ordered = resolveSampleOrder(corpus, beamtimeNum, loupeSampleId, sampleOrderState);
    const { index, prevId, nextId } = sampleNeighbors(ordered, loupeSampleId);
    if (index < 0 || ordered.length <= 1) return null;
    const goto = (id: number): void => {
      const params = new URLSearchParams();
      if (beamtimeNum !== undefined) params.set("beamtime", String(beamtimeNum));
      const qs = params.toString();
      navigate(`/samples/loupe/${id}${qs ? `?${qs}` : ""}`, { state: { sampleOrder: ordered } });
    };
    return (
      <SampleStepper
        name={sampleDisplayName(sample)}
        index={index}
        total={ordered.length}
        {...(prevId !== undefined ? { onPrev: () => goto(prevId) } : {})}
        {...(nextId !== undefined ? { onNext: () => goto(nextId) } : {})}
        prevKey={prevKey}
        nextKey={nextKey}
      />
    );
  })();

  // Compose from the print TopBar slot-shell (wordmark → children → spacer →
  // rightSlot). The shell owns the bar appearance (h-14, hairline, paper, the
  // 24px gutters, the flex-1 spacer); the corpus-topbar contract testid is
  // carried via TopBar's data-testid passthrough.
  return (
    <TopBar
      data-testid="corpus-topbar"
      wordmark={wordmark}
      rightSlot={focusStepper ?? loupeStepper}
    >
      {stageTabs}
      {beamtimeChip}
    </TopBar>
  );
}
