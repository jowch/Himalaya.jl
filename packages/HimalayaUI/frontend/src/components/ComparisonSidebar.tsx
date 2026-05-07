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
import { useComparisons } from "../queries";
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

  const sorted = useMemo(() => {
    return rows.slice().sort((a, b) => {
      // Most-recent first by `updated_at`. Null sorts last.
      const ai = a.updated_at ?? "";
      const bi = b.updated_at ?? "";
      if (ai === bi) return b.id - a.id; // tie-break by id desc
      if (ai === "") return 1;
      if (bi === "") return -1;
      return bi.localeCompare(ai);
    });
  }, [rows]);

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
            className="text-fg-muted text-sm p-3 text-center"
          >
            {rows.length === 0
              ? "No comparisons yet. Click + New to create one."
              : "No matches."}
          </li>
        ) : (
          filtered.map((c) => {
            const active = c.id === activeComparisonId;
            return (
              <li
                key={c.id}
                data-testid="comparison-list-item"
                data-comparison-id={c.id}
                {...(active ? { "data-active": "true" } : {})}
              >
                <button
                  type="button"
                  onClick={() => onPickRow(c.id)}
                  className={
                    "w-full text-left px-2 py-1.5 rounded text-sm "
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
              </li>
            );
          })
        )}
      </ul>
      </Skeleton>
    </div>
  );
}
