import { Link, useLocation, useSearchParams } from "react-router-dom";
import { useExperiments } from "../queries";

interface Stage {
  id: "samples" | "index" | "series";
  label: string;
  /** Absolute path the tab links to. Omitted = inert (disabled) tab. */
  to?: string;
}

// Samples (#160) and Series (#173) are live surfaces. Index stays inert until
// Phase 4 builds the focus workspace under a corpus path (redesign master
// plan §2.4). A tab with a `to` renders as a Link and derives its active state
// from the current route; a tab without one renders as a disabled button.
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
 * mockup omits them. The `T` theme shortcut (useGlobalShortcuts) still toggles
 * theme; the `username` identity state + persist survive untouched. A
 * multiplayer switch-user control + the theme-toggle button are deferred to a
 * follow-up issue (UI relocation, not in I5.1 scope).
 */
export function CorpusTopbar(): JSX.Element {
  const [searchParams, setSearchParams] = useSearchParams();
  const { pathname } = useLocation();
  const experimentsQuery = useExperiments();
  const beamtime = searchParams.get("beamtime") ?? "";

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
      <span
        data-testid="corpus-wordmark"
        className="text-sm font-bold uppercase tracking-[0.16em] text-ink"
      >
        Himalaya <span className="font-semibold text-ink-faint">· SAXS</span>
      </span>

      <nav data-testid="stage-tabs" aria-label="Workflow stages" className="flex gap-0.5">
        {STAGES.map((s) => {
          const dot = (
            <span
              aria-hidden="true"
              className={
                "inline-block w-1 h-1 rounded-full mr-1.5 align-middle " +
                (s.to !== undefined ? "bg-print-accent" : "bg-hair-strong")
              }
            />
          );
          // Active = this tab's path is the current route's prefix. Derived
          // from the router (not hardcoded) now that multiple stages are live.
          const isActive = s.to !== undefined && pathname.startsWith(s.to);
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
              disabled
              data-testid={`stage-tab-${s.id}`}
              className="px-2.5 py-1.5 rounded text-xs font-semibold uppercase
                         tracking-wide text-ink-faint cursor-not-allowed"
            >
              {dot}
              {s.label}
            </button>
          );
        })}
      </nav>

      <select
        data-testid="beamtime-chip"
        aria-label="Filter to a beamtime"
        value={beamtime}
        onChange={handlePick}
        className="rounded-full border border-hair-strong bg-plate px-2.5 py-1
                   text-xs font-semibold text-ink"
      >
        <option value="">Beamtime — all experiments</option>
        {(experimentsQuery.data ?? []).map((exp) => (
          <option key={exp.id} value={exp.id}>
            {exp.name ?? `Experiment ${exp.id}`}
          </option>
        ))}
      </select>

      <span className="flex-1" />
    </header>
  );
}
