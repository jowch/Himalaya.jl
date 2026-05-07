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
import type { ComparisonSummary } from "../api";
import { HintText } from "./ui";

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
  // recently pinned at top), then non-pinned by updated_at desc.
  const sorted = useMemo(() => {
    const pinnedOrder = pinsQ.data ?? [];
    const pinnedIdx = new Map<number, number>(
      pinnedOrder.map((id, i) => [id, i] as const),
    );
    const byUpdated = (a: ComparisonSummary, b: ComparisonSummary): number => {
      const ai = a.updated_at ?? "";
      const bi = b.updated_at ?? "";
      if (ai === bi) return b.id - a.id;
      if (ai === "") return 1;
      if (bi === "") return -1;
      return bi.localeCompare(ai);
    };
    const pinned: ComparisonSummary[] = [];
    const unpinned: ComparisonSummary[] = [];
    for (const c of rows) (pinSet.has(c.id) ? pinned : unpinned).push(c);
    pinned.sort((a, b) =>
      (pinnedIdx.get(a.id) ?? 0) - (pinnedIdx.get(b.id) ?? 0));
    unpinned.sort(byUpdated);
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
    if (experimentId !== undefined) navigate(`/experiments/${experimentId}/compare`);
  };
  const onPickAll = (): void => navigate("/compare/all");
  const onNew = (): void => {
    if (experimentId !== undefined) navigate(`/experiments/${experimentId}/compare/new`);
    else navigate("/compare/all"); // no experiment context → go to global list (new requires :eid)
  };
  const onPickRow = (id: number): void => {
    if (scope === "experiment" && experimentId !== undefined) {
      navigate(`/experiments/${experimentId}/compare/${id}`);
    } else {
      // Global scope → still navigate to the experiment-scoped review URL
      // is impossible without :eid; for now, stay on /compare/all and rely
      // on Phase 5's picker modal to land users on the right URL. We keep
      // the click as a no-op selector to surface the row visually.
      navigate(`/compare/all`);
    }
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
          disabled={experimentId === undefined}
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
      <ul className="flex-1 min-h-0 overflow-y-auto flex flex-col gap-1">
        {filtered.length === 0 ? (
          <li
            data-testid="comparison-sidebar-empty"
            className="text-fg-muted text-sm p-4 text-center flex flex-col items-center gap-2"
          >
            {rows.length === 0 ? (
              <>
                <span className="italic">No comparisons yet.</span>
                <button
                  type="button"
                  data-testid="comparison-sidebar-empty-new"
                  onClick={onNew}
                  disabled={experimentId === undefined}
                  className="px-3 py-1 rounded border border-border text-fg
                             hover:bg-bg-elevated text-sm disabled:opacity-50
                             disabled:cursor-not-allowed"
                >
                  + New comparison
                </button>
              </>
            ) : (
              <span>No matches.</span>
            )}
          </li>
        ) : (
          filtered.map((c) => {
            const active = c.id === activeComparisonId;
            const pinned = pinSet.has(c.id);
            const onTogglePin = (e: React.MouseEvent): void => {
              e.stopPropagation();
              if (pinned) unpin.mutate(c.id);
              else        pin.mutate(c.id);
            };
            return (
              <li
                key={c.id}
                data-testid="comparison-list-item"
                data-comparison-id={c.id}
                {...(active ? { "data-active": "true" } : {})}
                {...(pinned ? { "data-pinned": "true" } : {})}
                className="group flex items-start gap-1"
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
                  <div className="font-medium truncate">{c.title || `Comparison #${c.id}`}</div>
                  {c.description && (
                    <div className="text-xs text-fg-dim truncate">{c.description}</div>
                  )}
                </button>
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
          })
        )}
      </ul>
      </Skeleton>
    </div>
  );
}
