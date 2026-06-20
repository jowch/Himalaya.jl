import { Link, Outlet, useLocation, useParams } from "react-router-dom";
import { useState } from "react";
import { useExperiment } from "../../queries";
import { useAppState } from "../../state";
import * as api from "../../api";
import { authOpts } from "../../lib/authOpts";
import { Wordmark } from "../ui/Wordmark";
import { Kicker } from "../ui/Kicker";
import { Input } from "../ui/Input";
import { StatBar } from "../ui/StatBar";
import type { StatBarStat } from "../ui/StatBar";
import { ProgressBar } from "../ui/ProgressBar";
import { PageFrame } from "../components/PageFrame";
import { ExperimentTopNav } from "./ExperimentTopNav";

interface TabDef {
  id: "corpus" | "config";
  label: string;
}
const TABS: readonly TabDef[] = [
  { id: "corpus", label: "Corpus" },
  { id: "config", label: "Configuration" },
];

/**
 * ExperimentShell — the /experiments/:id layout route. Renders its OWN chrome
 * (wordmark + ExperimentTopNav + experiment header + Corpus|Configuration tab
 * bar) outside CorpusShell so the header never stacks on CorpusTopbar (spec
 * §7/§9.6). Children (ExperimentCorpusPage / ConfigurationPage) render in the
 * Outlet. The grouping-review route reuses this shell too but hides the tab
 * bar's active state (E2 wires the banner → grouping surface).
 */
export function ExperimentShell(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const expId = id ? Number(id) : 0;
  const { pathname } = useLocation();
  const exp = useExperiment(expId);
  const inFlight = useAppState((s) => s.ingestInFlight?.[expId]);
  const username = useAppState((s) => s.username);

  // Controlled draft for the edit-in-place name field. Initialized from the
  // loaded experiment; resets when the server data changes (e.g. after a
  // commit from another tab). The commit handler calls updateExperiment PATCH
  // (Task 1c widens the route to accept name as a plain-write field).
  const serverName = exp.data?.name ?? "";
  // Controlled draft for the edit-in-place name field.
  // `localEdit` tracks whether the user has typed locally since the last commit.
  // When false (no pending local edit), we render the current server name
  // directly so the Input always reflects the latest loaded value without a
  // separate effect/re-render cycle. When true, we render `pendingDraft`.
  const [localEdit, setLocalEdit] = useState(false);
  const [pendingDraft, setPendingDraft] = useState("");
  // The effective draft: if the user has typed, use their pending draft;
  // otherwise fall back to the server value (or the placeholder).
  const nameDraft = localEdit ? pendingDraft : serverName;

  const commitName = (): void => {
    if (!localEdit) return; // nothing pending
    const trimmed = pendingDraft.trim();
    setLocalEdit(false);
    setPendingDraft("");
    if (trimmed === serverName) return; // no-op
    // authOpts omits undefined keys (exactOptionalPropertyTypes: true).
    void api.updateExperiment(expId, { name: trimmed }, authOpts(username, undefined));
  };

  const status = inFlight?.status ?? exp.data?.ingest_status ?? "idle";
  const isProcessing = status === "scanning" || status === "analyzing";

  const stats: StatBarStat[] = isProcessing
    ? [
        { key: "processed", caption: "Processed",
          value: inFlight ? `${inFlight.processed} / ~${inFlight.total}` : "—" },
        { key: "span", caption: "Span", value: "pending", pending: true },
      ]
    : [
        // E1 ledger reads the experiment detail. Real counts (loads/samples/
        // exposures/sessions) arrive from the extended GET /api/experiments/:id
        // roll-up (spec §9.2); until then show placeholders the rollup fills.
        { key: "loads", caption: "Loads", value: "—" },
        { key: "samples", caption: "Samples", value: "—" },
        { key: "exposures", caption: "Exposures", value: "—" },
        { key: "sessions", caption: "Sessions", value: "—" },
      ];

  return (
    <div
      data-testid="experiment-shell"
      className="h-full w-full flex flex-col min-h-0 bg-paper text-ink overflow-auto"
    >
      {/* Own top chrome (NOT CorpusTopbar). */}
      <header className="flex items-center gap-6 px-6 h-14 border-b border-hair shrink-0">
        <Link
          to="/experiments"
          className="rounded focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent focus-visible:outline-offset-2"
        >
          <Wordmark tail="SAXS">Himalaya</Wordmark>
        </Link>
        <ExperimentTopNav />
      </header>

      <PageFrame width="experiment" className="px-6 py-6 flex-1 min-h-0 flex flex-col">
        {/* Experiment header */}
        <div className="flex items-start justify-between gap-6">
          <div className="min-w-0">
            <Kicker>Experiment</Kicker>
            {/* Edit-in-place name: Input variant='title' controlled by a local
                draft; commits on blur / Enter via updateExperiment PATCH (Task
                1c widens the route to accept name as a plain-write field).
                The wrapper is suppressed during loading so `findByTestId`
                waits for data to arrive before resolving. */}
            {exp.data !== undefined && (
              <Input
                variant="title"
                testId="experiment-header-name"
                value={nameDraft || `Experiment ${expId}`}
                onValueChange={(v) => {
                  setLocalEdit(true);
                  setPendingDraft(v);
                }}
                onBlur={commitName}
                onKeyDown={(e) => {
                  if (e.key === "Enter") {
                    commitName();
                    (e.target as HTMLElement).blur();
                  }
                }}
                aria-label="Experiment name"
              />
            )}
            {exp.data?.data_dir && (
              <p className="text-sm text-ink-soft font-mono mt-1 truncate">
                {exp.data.data_dir}
              </p>
            )}
          </div>
          <div
            data-testid="experiment-rescan-status"
            className="text-sm text-ink-soft shrink-0"
          >
            {isProcessing
              ? status === "scanning"
                ? "Scanning…"
                : "Analyzing…"
              : status === "failed"
                ? "Scan failed"
                : "Up to date"}
          </div>
        </div>

        {isProcessing && (
          <div className="mt-3">
            <ProgressBar
              value={inFlight ? inFlight.processed : 0}
              total={inFlight ? Math.max(inFlight.total, 1) : 1}
              label="Ingest progress"
            />
          </div>
        )}

        <StatBar aria-label="Experiment stats" stats={stats} className="my-5" />

        {/* Corpus | Configuration tab bar */}
        <nav
          data-testid="experiment-tab-bar"
          aria-label="Experiment views"
          className="flex gap-1 border-b border-hair"
        >
          {TABS.map((t) => {
            const to = `/experiments/${expId}/${t.id}`;
            const isActive = pathname.startsWith(to);
            return (
              <Link
                key={t.id}
                to={to}
                data-testid={`exp-tab-${t.id}`}
                aria-current={isActive ? "page" : undefined}
                className={
                  "px-3 py-2 text-sm font-semibold no-underline -mb-px border-b-2 " +
                  (isActive ? "text-ink border-accent" : "text-ink-soft border-transparent")
                }
              >
                {t.label}
              </Link>
            );
          })}
        </nav>

        <div className="flex-1 min-h-0 pt-5">
          <Outlet />
        </div>
      </PageFrame>
    </div>
  );
}
