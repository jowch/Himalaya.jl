import { Link, useSearchParams } from "react-router-dom";
import { useExperiments } from "../queries";

interface Stage {
  id: "samples" | "index" | "series";
  label: string;
  /** Absolute path the tab links to. Omitted = inert (disabled) tab. */
  to?: string;
}

// During the Phase-1 migration only the Samples surface exists in the new
// shell. Index and Series are inert tabs until Phases 4 and 3 build those
// surfaces (redesign master plan §2.4).
const STAGES: readonly Stage[] = [
  { id: "samples", label: "Samples", to: "/samples" },
  { id: "index", label: "Index" },
  { id: "series", label: "Series" },
];

/**
 * CorpusTopbar — the topbar for the redesigned corpus-scoped shell: a
 * corpus-level wordmark, the three workflow stage-tabs, and the Beamtime
 * facet chip.
 *
 * The Beamtime chip is an experiment picker: it reads and writes the
 * `?beamtime=<experiment_id>` URL query that the /samples contact sheet
 * (#160) filters on. The URL is the only channel — no prop coupling.
 */
export function CorpusTopbar(): JSX.Element {
  const [searchParams, setSearchParams] = useSearchParams();
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
          return s.to !== undefined ? (
            <Link
              key={s.id}
              to={s.to}
              data-testid={`stage-tab-${s.id}`}
              // Always active in #155 — Samples is the only live stage-tab.
              // When Index/Series go live (#3.x/#4.x), derive this from the route.
              data-active="true"
              className="px-2.5 py-1.5 rounded text-xs font-semibold uppercase
                         tracking-wide text-ink bg-paper-sunk no-underline"
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
