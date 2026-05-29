import { Link, useLocation, useNavigate, useSearchParams } from "react-router-dom";
import { useAppState } from "../state";
import { useCorpusSamples, useExperiments } from "../queries";
import { sampleDisplayName } from "../lib/sample/displayName";
import { Kicker } from "./ui/Kicker";

interface Stage {
  id: "samples" | "index" | "series";
  label: string;
  /** Absolute path the tab links to. Omitted = inert (disabled) tab. */
  to?: string;
}

// Samples (#160) and Series (#173) are live surfaces. The Index stage is the
// focus workspace at `/sample/:id` (Phase 4 cutover, #181) — but there is no
// canonical `/index` landing route (those URLs redirect), so the tab is not a
// Link. It renders as a button that is INERT off a sample route and ACTIVE
// (not navigable) on `/sample/:id` (R5 / #228, F-14). A tab with a `to`
// renders as a Link and derives its active state from the current route.
const STAGES: readonly Stage[] = [
  { id: "samples", label: "Samples", to: "/samples" },
  { id: "index", label: "Index" },
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
  // stepper, the active Index tab, and the Notes toggle only appear there.
  const onSampleRoute = pathname.startsWith("/sample/");
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const toggleNotesDrawer = useAppState((s) => s.toggleNotesDrawer);
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

  // F-12: the Notes toggle is the < xl fallback for the focus Notes margin. The
  // badge reflects whether the active sample carries notes (so the user knows
  // there's something behind the drawer without opening it).
  const hasNotes = (activeSample?.notes ?? "").trim().length > 0;

  function handlePick(event: React.ChangeEvent<HTMLSelectElement>): void {
    const value = event.target.value;
    setSearchParams((prev) => {
      const next = new URLSearchParams(prev);
      if (value === "") next.delete("beamtime");
      else next.set("beamtime", value);
      return next;
    });
  }

  return (
    <header
      data-testid="corpus-topbar"
      className="h-14 shrink-0 flex items-center gap-4 px-6 border-b border-hair bg-paper"
    >
      <Link
        to="/samples"
        data-testid="corpus-wordmark"
        aria-label="Himalaya SAXS, go to the corpus"
        className="text-sm font-bold uppercase tracking-[0.16em] text-ink rounded focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent focus-visible:outline-offset-2"
      >
        Himalaya <span className="font-semibold text-ink-faint">· SAXS</span>
      </Link>

      <nav data-testid="stage-tabs" aria-label="Workflow stages" className="flex gap-0.5">
        {STAGES.map((s) => {
          // The Index tab has no `to` but is live on a sample route (F-14), so
          // its dot lights when active too — not just for `to`-bearing tabs.
          const live = s.to !== undefined || (s.id === "index" && onSampleRoute);
          const dot = (
            <span
              aria-hidden="true"
              className={
                "inline-block w-1 h-1 rounded-full mr-1.5 align-middle " +
                (live ? "bg-print-accent" : "bg-hair-strong")
              }
            />
          );
          // Active = this tab's path is the current route's prefix. Derived
          // from the router (not hardcoded) now that multiple stages are live.
          // F-14: the linkless Index tab is active on `/sample/:id` (the focus
          // workspace), reflecting the focus surface as the Index stage.
          const isActive = s.to !== undefined
            ? pathname.startsWith(s.to)
            : s.id === "index" && onSampleRoute;
          return s.to !== undefined ? (
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
              {dot}
              {s.label}
            </Link>
          ) : (
            <button
              key={s.id}
              type="button"
              // Active (on /sample/:id) the Index tab is a non-navigable
              // current-stage marker, not a disabled control. Off a sample
              // route it stays inert.
              disabled={!isActive}
              data-testid={`stage-tab-${s.id}`}
              data-active={isActive ? "true" : undefined}
              aria-current={isActive ? "page" : undefined}
              className={
                "px-2.5 py-1.5 rounded text-xs font-semibold uppercase " +
                "tracking-wide " +
                (isActive
                  ? "text-ink bg-paper-sunk cursor-default"
                  : "text-ink-faint cursor-not-allowed")
              }
            >
              {dot}
              {s.label}
            </button>
          );
        })}
      </nav>

      {/* The beamtime filter is honored only on the samples surface (it filters
          the contact sheet, and the loupe back-link preserves it). On /series
          and /sample/:id it was changeable but inert — hide it there rather
          than present a control that does nothing. */}
      {pathname.startsWith("/samples") && (
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
      )}

      <span className="flex-1" />

      {/* M-9: Contact sheet | Loupe segmented view switch. The loupe is a
          sample-scoped route (/samples/loupe/:id); the topbar has no sample id
          on /samples, so the switch is route-based — on a loupe route both
          segments reflect state and "Loupe" is active; on the sheet, "Contact
          sheet" is active and "Loupe" is disabled (no sample is selected to
          open). The Contact-sheet link preserves the ?beamtime= filter.

          R5 (#228) coordination: scoped to the samples surface. The focus
          surface (/sample/:id) shows its own per-sample stepper + Notes toggle
          below, not the contact-sheet/loupe switch (which is meaningless
          there). */}
      {pathname.startsWith("/samples") && (() => {
        const onLoupe = pathname.startsWith("/samples/loupe");
        const beamtimeQuery = beamtime === "" ? "" : `?beamtime=${beamtime}`;
        const sheetHref = `/samples${beamtimeQuery}`;
        const segBase =
          "px-3 py-1.5 text-[11.5px] font-semibold no-underline";
        const active = "bg-ink text-paper";
        const inactive = "text-ink-faint";
        return (
          <div
            data-testid="view-seg"
            className="flex overflow-hidden rounded-md border border-hair-strong bg-plate"
          >
            <Link
              to={sheetHref}
              data-testid="view-seg-sheet"
              data-active={onLoupe ? undefined : "true"}
              aria-current={onLoupe ? undefined : "page"}
              className={`${segBase} ${onLoupe ? inactive : active}`}
            >
              Contact sheet
            </Link>
            {onLoupe ? (
              <span
                data-testid="view-seg-loupe"
                data-active="true"
                aria-current="page"
                className={`${segBase} ${active}`}
              >
                Loupe
              </span>
            ) : (
              <button
                type="button"
                disabled
                data-testid="view-seg-loupe"
                title="Open a sample to use the loupe"
                className={`${segBase} ${inactive} cursor-not-allowed bg-transparent`}
              >
                Loupe
              </button>
            )}
          </div>
        );
      })()}

      {/* F-13: per-sample stepper — the focus surface's primary inter-sample
          nav. URL-routed so the one-way route→store sync stays intact. */}
      {showStepper && (
        <div
          data-testid="sample-stepper"
          className="flex items-center gap-2 text-ink"
        >
          <button
            type="button"
            data-testid="sample-stepper-prev"
            disabled={prevSample === undefined}
            onClick={() => prevSample && navigate(`/sample/${prevSample.id}`)}
            aria-label="Previous sample"
            className="rounded px-1.5 py-0.5 text-base leading-none text-ink-faint
                       hover:text-ink disabled:cursor-not-allowed disabled:opacity-30"
          >
            &#8249;
          </button>
          <span className="flex flex-col items-end leading-tight">
            <span className="text-xs font-semibold text-ink">
              {sampleDisplayName(activeSample!)}
            </span>
            <Kicker as="span" tone="faint">sample {stepIdx + 1} of {siblings.length}</Kicker>
          </span>
          <button
            type="button"
            data-testid="sample-stepper-next"
            disabled={nextSample === undefined}
            onClick={() => nextSample && navigate(`/sample/${nextSample.id}`)}
            aria-label="Next sample"
            className="rounded px-1.5 py-0.5 text-base leading-none text-ink-faint
                       hover:text-ink disabled:cursor-not-allowed disabled:opacity-30"
          >
            &#8250;
          </button>
        </div>
      )}

      {/* F-12: Notes toggle — the < xl fallback for the focus Notes margin.
          Hidden at xl+ (mockup `#notes-btn{display:none}` above 1320px) where
          the persistent Notes margin column is shown instead; below xl the
          margin yields to this button + the drawer it toggles. */}
      {onSampleRoute && (
        <button
          type="button"
          data-testid="notes-toggle"
          data-has-notes={hasNotes ? "true" : undefined}
          onClick={toggleNotesDrawer}
          className="xl:hidden flex items-center gap-1.5 rounded border border-hair-strong
                     bg-plate px-2.5 py-1 text-xs font-semibold text-ink
                     hover:border-ink-faint"
        >
          Notes
          {hasNotes && (
            <span
              aria-hidden="true"
              className="inline-block h-1.5 w-1.5 rounded-full bg-print-accent"
            />
          )}
        </button>
      )}
    </header>
  );
}
