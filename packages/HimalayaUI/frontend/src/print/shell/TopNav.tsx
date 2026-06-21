import { Link, useLocation } from "react-router-dom";
import { TopBar } from "../ui/TopBar";
import { Wordmark } from "../ui/Wordmark";

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
 * TopNav — the unified top navigation bar (App Shell Unification spec §3.1).
 * Wordmark "Himalaya · SAXS" links to /experiments; two section tabs
 * (Experiments, Series) with router-derived active state. Replaces both
 * CorpusTopbar and ExperimentTopNav in T3.2 (deletion step). No Samples tab,
 * no Beamtime chip, no gear (Configuration rides the experiment header, §3.1).
 */
export function TopNav(): JSX.Element {
  const { pathname } = useLocation();

  const wordmark = (
    <Link
      to="/experiments"
      data-testid="topnav-wordmark"
      aria-label="Himalaya SAXS, go to experiments"
      className="rounded focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent focus-visible:outline-offset-2"
    >
      <Wordmark tail="SAXS">Himalaya</Wordmark>
    </Link>
  );

  const tabs = (
    <nav data-testid="topnav-tabs" aria-label="Sections" className="flex gap-0.5">
      {ITEMS.map((it) => {
        const isActive = pathname.startsWith(it.to);
        return (
          <Link
            key={it.id}
            to={it.to}
            data-testid={`topnav-tab-${it.id}`}
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

  return (
    <TopBar data-testid="topnav" wordmark={wordmark}>
      {tabs}
    </TopBar>
  );
}
