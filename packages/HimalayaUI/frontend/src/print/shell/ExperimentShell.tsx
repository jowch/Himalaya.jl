import { Outlet, useNavigate, useParams } from "react-router-dom";
import { useState } from "react";
import { useExperiment, useTriggerScan } from "../../queries";
import { useAppState } from "../../state";
import * as api from "../../api";
import { authOpts } from "../../lib/authOpts";
import { Kicker } from "../ui/Kicker";
import { Input } from "../ui/Input";
import { Button } from "../ui/Button";
import { IconButton } from "../ui/IconButton";
import { StatBar } from "../ui/StatBar";
import type { StatBarStat } from "../ui/StatBar";
import { ProgressBar } from "../ui/ProgressBar";
import { PageFrame } from "../components/PageFrame";
import { effectiveIngestStatus } from "../../lib/ingestStatus";

/**
 * ExperimentShell — the /experiments/:id layout route. Renders the experiment
 * page content (header + Outlet). T3.2: the top chrome (TopNav) is provided by
 * the outer AppShell — this is pure page content, not a chrome. The
 * Corpus|Configuration tab bar is retired (M3): Corpus is the experiment home
 * (the index route), and Configuration is reached via the ⚙ in the header.
 */
export function ExperimentShell(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const expId = id ? Number(id) : 0;
  const navigate = useNavigate();
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

  // Terminal persisted status overrides a stale "scanning" overlay (8c) — see
  // effectiveIngestStatus; useExperiment self-heals the persisted value.
  const status = effectiveIngestStatus(inFlight?.status, exp.data?.ingest_status);
  const isProcessing = status === "scanning" || status === "analyzing";

  const triggerScan = useTriggerScan(expId);

  const expStats = exp.data?.stats;
  const stats: StatBarStat[] = isProcessing
    ? [
        { key: "processed", caption: "Processed",
          value: inFlight ? `${inFlight.processed} / ~${inFlight.total}` : "—" },
        { key: "span", caption: "Span", value: "pending", pending: true },
      ]
    : [
        // E1 ledger reads the experiment detail. Real counts arrive from the
        // stats roll-up on GET /api/experiments/:id (spec §9.2).
        { key: "loads", caption: "Loads", value: expStats != null ? String(expStats.loads) : "—" },
        { key: "samples", caption: "Samples", value: expStats != null ? String(expStats.samples) : "—" },
        { key: "exposures", caption: "Exposures", value: expStats != null ? String(expStats.exposures) : "—" },
        { key: "span", caption: "Span", value: expStats != null ? `${expStats.span_hours}h` : "—" },
        { key: "sessions", caption: "Sessions", value: expStats != null ? String(expStats.sessions) : "—" },
      ];

  return (
    <div
      data-testid="experiment-shell"
      className="h-full w-full flex flex-col min-h-0 bg-paper text-ink overflow-auto"
    >
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
          <div className="flex items-center gap-2 shrink-0">
            <div
              data-testid="experiment-rescan-status"
              className="text-sm text-ink-soft"
            >
              {isProcessing
                ? status === "scanning"
                  ? "Scanning…"
                  : "Analyzing…"
                : status === "failed"
                  ? "Scan failed"
                  : "Up to date"}
            </div>
            <Button
              variant="outline"
              data-testid="experiment-rescan-button"
              disabled={isProcessing}
              onClick={() => triggerScan.mutate()}
            >
              Rescan
            </Button>
            <IconButton
              label="Configuration"
              tone="ghost"
              data-testid="experiment-config-gear"
              onClick={() => navigate(`/experiments/${expId}/config`)}
            >
              ⚙
            </IconButton>
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

        <StatBar aria-label="Experiment stats" stats={stats} className="mt-5 mb-2" />

        <div className="flex-1 min-h-0 pt-5">
          <Outlet />
        </div>
      </PageFrame>
    </div>
  );
}
