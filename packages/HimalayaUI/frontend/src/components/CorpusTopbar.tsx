import { Link } from "react-router-dom";

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
 * corpus-level wordmark (experiment is a facet now, not the crumb), the
 * three workflow stage-tabs, and the Beamtime facet chip.
 *
 * The Beamtime chip is a presentational placeholder in #155 — `?beamtime=`
 * URL query state is owned by the `/samples` route (#160).
 */
export function CorpusTopbar(): JSX.Element {
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

      <nav data-testid="stage-tabs" className="flex gap-0.5">
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

      <button
        type="button"
        data-testid="beamtime-chip"
        title="Filter to a beamtime (coming soon)"
        className="flex items-center gap-1.5 rounded-full border border-hair-strong
                   bg-plate px-2.5 py-1 text-xs font-semibold text-ink"
      >
        <span>Beamtime</span>
        <span aria-hidden="true" className="text-ink-faint">
          ▾
        </span>
      </button>

      <span className="flex-1" />
    </header>
  );
}
