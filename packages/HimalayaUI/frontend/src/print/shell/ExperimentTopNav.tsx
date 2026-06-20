import { Link, useLocation } from "react-router-dom";

interface NavItem {
  id: "experiments" | "series";
  label: string;
  to: string;
}

const ITEMS: readonly NavItem[] = [
  { id: "experiments", label: "Experiments", to: "/experiments" },
  { id: "series", label: "Series", to: "/series" },
];

/**
 * ExperimentTopNav — the two top-level axes (spec §7: Experiments | Series).
 * Used by ExperimentShell, which renders its own chrome OUTSIDE CorpusShell so
 * the experiment header never stacks on CorpusTopbar. Router Links with
 * `startsWith` active detection (mirrors CorpusTopbar's stage tabs). The home
 * phase bar is intentionally dropped (round-3 note).
 */
export function ExperimentTopNav(): JSX.Element {
  const { pathname } = useLocation();
  return (
    <nav data-testid="experiment-top-nav" aria-label="Sections" className="flex gap-0.5">
      {ITEMS.map((it) => {
        const isActive = pathname.startsWith(it.to);
        return (
          <Link
            key={it.id}
            to={it.to}
            data-testid={`topnav-${it.id}`}
            aria-current={isActive ? "page" : undefined}
            className={
              "px-2.5 py-1.5 rounded text-xs font-semibold uppercase tracking-wide no-underline " +
              (isActive ? "text-ink bg-paper-sunk" : "text-ink-soft")
            }
          >
            <span
              aria-hidden="true"
              className="inline-block w-1 h-1 rounded-full mr-1.5 align-middle bg-print-accent"
            />
            {it.label}
          </Link>
        );
      })}
    </nav>
  );
}
