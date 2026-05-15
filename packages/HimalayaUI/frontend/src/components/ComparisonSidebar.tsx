/**
 * ComparisonSidebar (Plan §Phase 4, Task 4.2).
 *
 * Lists comparisons for the active scope (this experiment / all
 * experiments), sorted most-recent first. Hosts the scope toggle, search
 * box, and "+ New" button. Active comparison is marked via
 * `data-active="true"` so the page CSS can highlight it without component
 * coupling.
 *
 * Scope toggle behaviour:
 *   - "This experiment": stays at `/experiments/:eid/compare(/:id)`
 *   - "All experiments": navigates to `/compare/all`
 * Pin toggle is deferred to Phase 13 per the plan.
 */
import { useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import {
  useComparisons, useComparisonPins, usePinComparison, useUnpinComparison,
} from "../queries";
import { comparePath } from "../lib/comparison/routes";
import { relativeTime } from "../lib/comparison/relativeTime";
import { useCurrentUserId } from "../hooks/useCurrentUserId";
import type { ComparisonSummary } from "../api";
import { HintText } from "./ui";
import { useAppState } from "../state";

/**
 * Phase-summary line for a sidebar row (Compare UX F-1). Shows up to three
 * distinct member phases, an "+N more" overflow, and the trace count:
 *   ["Pn3m","Hex","Lam"], 4 → "Pn3m · Hex · Lam · 4 traces"
 *   five phases,          5 → "Pn3m · Im3m · Ia3d · +2 more · 5 traces"
 *   [],                   2 → "2 traces"
 */
function formatPhaseSummary(phases: string[], total: number): string {
  if (phases.length === 0) return `${total} traces`;
  if (phases.length <= 3) return `${phases.join(" · ")} · ${total} traces`;
  return `${phases.slice(0, 3).join(" · ")} · +${phases.length - 3} more · ${total} traces`;
}

// Mock fixture for boneyard layout capture. Renders a few canonical rows so
// the captured bones reflect the realistic geometry the user will see during
// a true cold fetch.
const COMPARISON_SIDEBAR_FIXTURE = (
  <ul className="flex-1 min-h-0 overflow-y-auto flex flex-col gap-1">
    {[
      { id: 1, title: "Pn3m vs Im3m sweep" },
      { id: 2, title: "Hex run B comparison" },
      { id: 3, title: "Lipid:water ratio" },
    ].map((c) => (
      <li key={c.id}>
        <div className="w-full text-left px-2 py-1.5 rounded text-sm text-fg-muted">
          <div className="font-medium truncate">{c.title}</div>
          <div className="text-xs text-fg-dim truncate">3 traces</div>
        </div>
      </li>
    ))}
  </ul>
);

interface Props {
  experimentId: number | undefined;
  scope: "experiment" | "all";
  activeComparisonId: number | undefined;
}

export function ComparisonSidebar({
  experimentId,
  scope,
  activeComparisonId,
}: Props): JSX.Element {
  const navigate = useNavigate();
  const [search, setSearch] = useState("");
  const currentUserId = useCurrentUserId();
  // Server id of the comparison currently open as an unsaved draft, if any.
  // Undefined for a brand-new (never-saved) draft — it has no listing row.
  const draftId = useAppState((s) => s.activeDraft?.id);

  // Fall back to the persisted Zustand experiment context when the URL has
  // no `:eid` (e.g. on `/compare/all`). Without this fallback the
  // "This experiment" tab is stuck-disabled the moment the user toggles
  // to All experiments — even though they have a live experiment context
  // they were just looking at (issue #79).
  const fallbackExperimentId = useAppState((s) => s.activeExperimentId);
  const effectiveExperimentId = experimentId ?? fallbackExperimentId;

  // The query key picks up the right scope. When scope is "experiment" but
  // experimentId is undefined (shouldn't happen in practice), fall back to
  // "all" so we never pass a bad id to the experiment-listing route.
  const listScope: number | "all" =
    scope === "all" || experimentId === undefined ? "all" : experimentId;
  const listQ = useComparisons(listScope);
  const rows = listQ.data ?? [];

  // Phase 13 — per-user pins. The pin set is a Set<number> for O(1) lookup
  // when partitioning the sorted list. Pinned ids that aren't in the
  // current scope's listing are silently ignored (a pinned comparison from
  // another experiment won't appear under "This experiment").
  const pinsQ = useComparisonPins();
  const pinSet = useMemo(
    () => new Set<number>(pinsQ.data ?? []),
    [pinsQ.data],
  );
  const pin   = usePinComparison();
  const unpin = useUnpinComparison();

  // Sort: pinned first (preserving the pin order from the API — most
  // recently pinned at top), then non-pinned by `last_event_at` desc
  // (spec §8.4 — replaces `updated_at`; `last_event_at` also covers chat
  // activity). ISO strings from SQLite are canonical `YYYY-MM-DDTHH:MM:SSZ`,
  // so lexicographic compare == chronological. Nulls sort LAST (a row with
  // no events is conceptually older than any timestamped row).
  const sorted = useMemo(() => {
    const pinnedOrder = pinsQ.data ?? [];
    const pinnedIdx = new Map<number, number>(
      pinnedOrder.map((id, i) => [id, i] as const),
    );
    const byLastEvent = (a: ComparisonSummary, b: ComparisonSummary): number => {
      const at = a.last_event_at;
      const bt = b.last_event_at;
      if (at === null && bt === null) return b.id - a.id;
      if (at === null) return 1;   // a is null → after b
      if (bt === null) return -1;  // b is null → after a
      if (at === bt) return b.id - a.id;
      return bt.localeCompare(at);
    };
    const pinned: ComparisonSummary[] = [];
    const unpinned: ComparisonSummary[] = [];
    for (const c of rows) (pinSet.has(c.id) ? pinned : unpinned).push(c);
    pinned.sort((a, b) =>
      (pinnedIdx.get(a.id) ?? 0) - (pinnedIdx.get(b.id) ?? 0));
    unpinned.sort(byLastEvent);
    return [...pinned, ...unpinned];
  }, [rows, pinSet, pinsQ.data]);

  const filtered = useMemo(() => {
    const q = search.trim().toLowerCase();
    if (q === "") return sorted;
    return sorted.filter((c: ComparisonSummary) => {
      const t = (c.title ?? "").toLowerCase();
      const d = (c.description ?? "").toLowerCase();
      return t.includes(q) || d.includes(q);
    });
  }, [sorted, search]);

  const onPickThis = (): void => {
    // Read from the URL prop OR the Zustand fallback — clicking "This
    // experiment" while on /compare/all should still navigate to the
    // user's last-active experiment context (issue #79).
    if (effectiveExperimentId !== undefined) {
      navigate(`/experiments/${effectiveExperimentId}/compare`);
    }
  };
  const onPickAll = (): void => navigate("/compare/all");
  const onNew = (): void => {
    // From the global scope we now have a real /compare/all/new route, so
    // creating without an experiment context lands somewhere coherent
    // (was previously a falsey navigate to /compare/all).
    navigate(
      comparePath({
        scope: scope === "experiment" && experimentId !== undefined
          ? "experiment"
          : "all",
        eid: experimentId,
        isNew: true,
      }),
    );
  };
  const onPickRow = (id: number): void => {
    // Both scopes deep-link directly to the comparison's review page now
    // — global scope routes to /compare/all/:id (Phase 4 follow-up); the
    // earlier "stay on /compare/all and let the picker land you" hack is gone.
    navigate(
      comparePath({
        scope: scope === "experiment" && experimentId !== undefined
          ? "experiment"
          : "all",
        eid: experimentId,
        id,
      }),
    );
  };

  const dataScope: "this" | "all" = scope === "experiment" ? "this" : "all";

  return (
    <div
      data-testid="comparison-sidebar"
      className="flex-1 min-h-0 flex flex-col gap-2 p-3"
    >
      <div
        data-testid="comparison-scope-toggle"
        data-scope={dataScope}
        role="tablist"
        className="inline-flex items-center gap-0.5 p-0.5
                   bg-bg-elevated border border-border rounded-full text-xs"
      >
        <button
          type="button"
          role="tab"
          data-testid="comparison-scope-this"
          aria-selected={dataScope === "this"}
          onClick={onPickThis}
          disabled={effectiveExperimentId === undefined}
          className={
            "px-2.5 py-0.5 rounded-full font-medium "
            + (dataScope === "this"
              ? "bg-accent/15 text-accent"
              : "text-fg-muted hover:text-fg")
          }
        >
          This experiment
        </button>
        <button
          type="button"
          role="tab"
          data-testid="comparison-scope-all"
          aria-selected={dataScope === "all"}
          onClick={onPickAll}
          className={
            "px-2.5 py-0.5 rounded-full font-medium "
            + (dataScope === "all"
              ? "bg-accent/15 text-accent"
              : "text-fg-muted hover:text-fg")
          }
        >
          All experiments
        </button>
      </div>

      <div className="flex items-center gap-2">
        <input
          data-testid="comparison-sidebar-search"
          type="search"
          placeholder="Search comparisons"
          value={search}
          onChange={(e) => setSearch(e.target.value)}
          className="flex-1 bg-transparent border border-border rounded px-2 py-1 text-sm"
        />
        <button
          type="button"
          data-testid="comparison-new"
          onClick={onNew}
          className="px-2 py-1 rounded border border-border text-fg-muted hover:text-fg text-sm"
        >
          + New
        </button>
      </div>

      <Skeleton
        name="comparison-sidebar"
        className="flex-1 min-h-0 flex flex-col"
        loading={listQ.isLoading}
        stagger={50}
        transition={200}
        fixture={COMPARISON_SIDEBAR_FIXTURE}
        fallback={<div className="p-3"><HintText>Loading comparisons…</HintText></div>}
      >
      {rows.length === 0 ? (
        <div
          data-testid="comparison-sidebar-empty"
          className="px-4 py-8 text-center text-fg-muted
                     flex flex-col items-center gap-3"
        >
          <p className="text-sm">
            {scope === "experiment"
              ? "No comparisons in this experiment yet."
              : "No comparisons yet."}
          </p>
          <button
            type="button"
            data-testid="sidebar-empty-new"
            onClick={onNew}
            className="px-3 py-1 rounded border border-border text-fg
                       hover:bg-bg-elevated text-sm"
          >
            + New comparison
          </button>
        </div>
      ) : filtered.length === 0 ? (
        <div
          data-testid="comparison-sidebar-no-matches"
          className="px-4 py-8 text-center text-fg-muted
                     flex flex-col items-center gap-3"
        >
          <p className="text-sm">No matches for &ldquo;{search}&rdquo;.</p>
          <button
            type="button"
            data-testid="sidebar-empty-clear"
            onClick={() => setSearch("")}
            className="px-3 py-1 rounded border border-border text-fg
                       hover:bg-bg-elevated text-sm"
          >
            Clear search
          </button>
        </div>
      ) : (
        <ul className="flex-1 min-h-0 overflow-y-auto flex flex-col gap-1">
          {filtered.map((c) => {
            const active = c.id === activeComparisonId;
            const pinned = pinSet.has(c.id);
            const onTogglePin = (e: React.MouseEvent): void => {
              e.stopPropagation();
              if (pinned) unpin.mutate(c.id);
              else        pin.mutate(c.id);
            };
            // Author byline: "by you" when the current user authored it,
            // else "by <username>", else "by —" when the author is unknown.
            const isMine = currentUserId !== undefined
              && c.created_by === currentUserId;
            const byline = isMine
              ? "by you"
              : c.author_username !== null ? `by ${c.author_username}` : "by —";
            const rel = relativeTime(c.last_event_at, Date.now());
            // This row has an unsaved draft open against it (F-2).
            const isDraft = draftId !== undefined && draftId === c.id;
            return (
              <li
                key={c.id}
                data-testid="comparison-list-item"
                data-comparison-id={c.id}
                {...(active ? { "data-active": "true" } : {})}
                {...(pinned ? { "data-pinned": "true" } : {})}
                className="group relative flex items-start gap-1"
              >
                <button
                  type="button"
                  onClick={() => onPickRow(c.id)}
                  className={
                    "flex-1 min-w-0 text-left px-2 py-1.5 rounded text-sm "
                    + (active
                      ? "bg-accent/10 text-fg"
                      : "text-fg-muted hover:text-fg hover:bg-bg-elevated")
                  }
                >
                  <div className="font-medium truncate">
                    {isDraft && (
                      <span
                        data-testid="sidebar-draft-dot"
                        aria-hidden="true"
                        className="text-accent mr-1"
                      >
                        •
                      </span>
                    )}
                    {(c.title || `Comparison #${c.id}`)
                      + (isDraft ? " (draft)" : "")}
                  </div>
                  <div className="text-xs text-fg-dim truncate">
                    {formatPhaseSummary(c.member_phases, c.member_count)}
                  </div>
                  <div className="text-xs text-fg-dim truncate">
                    {isDraft
                      ? "by you · just now"
                      : `${byline} · ${rel === null ? "—" : `edited ${rel}`}`}
                  </div>
                </button>
                {c.has_stale_members && (
                  <span
                    data-testid="sidebar-stale-warn"
                    title="Some members have stale indices"
                    className="absolute top-1 right-8 text-warning"
                  >
                    ⚠
                  </span>
                )}
                <button
                  type="button"
                  data-testid="comparison-pin-toggle"
                  data-comparison-id={c.id}
                  data-pinned={pinned ? "true" : "false"}
                  onClick={onTogglePin}
                  title={pinned ? "Unpin from top" : "Pin to top"}
                  aria-label={pinned ? "Unpin comparison" : "Pin comparison"}
                  className={
                    "shrink-0 px-1.5 py-1 rounded text-sm leading-none "
                    + (pinned
                      ? "text-accent"
                      : "text-fg-dim opacity-0 group-hover:opacity-100 hover:text-fg")
                  }
                >
                  {pinned ? "★" : "☆"}
                </button>
              </li>
            );
          })}
        </ul>
      )}
      </Skeleton>
    </div>
  );
}
