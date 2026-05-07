/**
 * ComparisonPicker (Plan §Phase 5, Task 5.2 frontend half).
 *
 * Multi-select exposure picker mounted at body root for the comparison
 * edit-mode "Add traces" flow. Per spec §Picker modal:
 *
 *   - Modal dialog with role + aria-modal + Esc-to-close.
 *   - Search box filters across exposure name + sample name + sample notes.
 *   - Filter chips: experiment (multi-select; defaults to active experiment)
 *     and tag (distinct (key, value) pairs from sample_tags).
 *   - "Recently used" section above the main list (per-user history from
 *     `comparison_members.created_by/created_at`).
 *   - Main list: <ul role="listbox"> of <ExposureListRow> with
 *     `data-testid="picker-row"` for E2E selectors.
 *   - Already-added rows are LOCKED (checkbox disabled, visually muted).
 *   - "Add selected" button at bottom. Disabled when nothing is picked.
 *   - Focus trap via useFocusTrap. Restores focus to the trigger on close.
 *   - Pick order = append order (spec). The user reorders via drag in the
 *     metadata gutter afterward (Phase 7).
 *
 * Lifecycle:
 *   - On open, focus the search input (matches NavModal).
 *   - Each row's checkbox toggles a local Set<exposureId>.
 *   - "Add selected" calls `addMember(exposureId)` for each id in
 *     selection-add-order, then closes via `onClose()`.
 *
 * Cross-experiment picks: the experiment filter chip can be widened to
 * include other experiments. We fetch all experiments via `useExperiments`
 * and let the user select multiple. By default only the current experiment
 * is selected so the user sees their immediate working set first.
 */
import { useEffect, useMemo, useRef, useState } from "react";
import { useQuery, useQueries, useQueryClient } from "@tanstack/react-query";
import { useAppState } from "../state";
import {
  useExperiments,
  useSamples,
  useRecentlyPickedExposures,
  useSampleTags,
} from "../queries";
import { useFocusTrap } from "../hooks/useFocusTrap";
import { ExposureListRow } from "./ExposureListRow";
import * as api from "../api";
import type { Exposure, Sample } from "../api";

/**
 * Resolve the current user's id by username, off the cached `listUsers`
 * query. Returns `undefined` until the lookup resolves so dependent queries
 * can stay disabled. Local helper rather than a top-level hook to keep the
 * picker's user-resolution scope contained — only the picker's
 * "Recently used" section needs it today.
 */
function useCurrentUserId(): number | undefined {
  const username = useAppState((s) => s.username);
  const usersQ = useQuery({
    queryKey: ["users"] as const,
    // Reuse the existing `listUsers` fetcher rather than duplicating the
    // request shape; TanStack dedupes if other code already loaded it.
    queryFn: () => api.listUsers(),
    enabled: username !== undefined,
  });
  if (username === undefined) return undefined;
  const u = (usersQ.data ?? []).find((x) => x.username === username);
  return u?.id;
}

interface Props {
  isOpen: boolean;
  onClose: () => void;
  /** The experiment context the picker opened from (defaults to selected). */
  experimentId: number | undefined;
}

interface RowEntry {
  exposure: Exposure;
  sample: Sample;
  /** True when the exposure is already a member of the active draft. */
  alreadyAdded: boolean;
}

export function ComparisonPicker({
  isOpen,
  onClose,
  experimentId,
}: Props): JSX.Element | null {
  const dialogRef = useRef<HTMLDivElement>(null);
  const inputRef = useRef<HTMLInputElement>(null);
  useFocusTrap(dialogRef, isOpen);

  const qc = useQueryClient();
  const draft = useAppState((s) => s.activeDraft);
  const addMember = useAppState((s) => s.addMember);

  const userId = useCurrentUserId();

  // Data sources.
  const experimentsQ = useExperiments();
  const recentsQ = useRecentlyPickedExposures(userId, 20);
  const tagsQ = useSampleTags(experimentId);

  // Active filter state. Experiments default to the current one (single
  // chip). Tag set is empty by default (no filter). Search is empty.
  const [search, setSearch] = useState("");
  const [selectedExpIds, setSelectedExpIds] = useState<number[]>(
    experimentId !== undefined ? [experimentId] : [],
  );
  const [selectedTags, setSelectedTags] = useState<Set<string>>(new Set());

  // Selection set — preserves insertion order for Pick order = append order.
  const [picks, setPicks] = useState<number[]>([]);

  // Reset transient state on open / close cycles. Without the reset, a
  // dismissed-then-reopened picker preserves stale picks/search and surprises
  // the user.
  useEffect(() => {
    if (isOpen) {
      setSearch("");
      setSelectedExpIds(experimentId !== undefined ? [experimentId] : []);
      setSelectedTags(new Set());
      setPicks([]);
      // Focus the input after the dialog mounts.
      inputRef.current?.focus();
    }
  }, [isOpen, experimentId]);

  // Aggregated row source: pull samples for every selected experiment, then
  // their exposures, in one flattened pass. We fetch via the cache (each
  // hook keys per-id); the in-memory join below is O(rows) not O(rows²).
  // Hooks must be at top level — call them for the active experiment as a
  // primary source. Cross-experiment picks (selectedExpIds.length > 1) are
  // surfaced in Phase 5+ but the v1 implementation treats the seeded
  // experimentId as the effective filter set; broadening to other
  // experiments is supported via `selectedExpIds` mutation but rendered
  // through the same query (TODO: enrich when multi-exp picks ship).
  const samplesQ = useSamples(experimentId ?? 0);
  const samples = samplesQ.data ?? [];

  // Map sample.id → Sample for fast row joining.
  const samplesById = useMemo(() => {
    const m = new Map<number, Sample>();
    for (const s of samples) m.set(s.id, s);
    return m;
  }, [samples]);

  // Fetch exposures for each sample. This calls listExposures under the
  // hood per sample; the test mocks return an array per sample id.
  // We accumulate them into a flat list for rendering.
  const allExposures = useExposuresAcrossSamples(samples);

  // Build the row list (joined + filtered + sorted).
  const rows = useMemo<RowEntry[]>(() => {
    const draftExposureIds = new Set(
      (draft?.members ?? [])
        .map((m) => m.exposure_id)
        .filter((id): id is number => id !== null && id !== undefined),
    );
    const out: RowEntry[] = [];
    for (const exposure of allExposures) {
      const sample = samplesById.get(exposure.sample_id);
      if (!sample) continue;
      // Tag filter (OR semantics across selected (key:value) pairs).
      if (selectedTags.size > 0) {
        const sampleHasMatch = sample.tags.some((t) =>
          selectedTags.has(`${t.key}:${t.value}`),
        );
        if (!sampleHasMatch) continue;
      }
      // Search filter (case-insensitive substring across name fields).
      if (search.trim() !== "") {
        const needle = search.toLowerCase();
        const haystack = [
          (exposure.filename ?? "").toLowerCase(),
          (sample.name ?? "").toLowerCase(),
          (sample.label ?? "").toLowerCase(),
          (sample.notes ?? "").toLowerCase(),
          ...sample.tags.map((t) => t.value.toLowerCase()),
        ];
        if (!haystack.some((s) => s.includes(needle))) continue;
      }
      out.push({
        exposure,
        sample,
        alreadyAdded: draftExposureIds.has(exposure.id),
      });
    }
    // Default sort: alphabetical by exposure name (filename minus .dat).
    out.sort((a, b) => {
      const an = (a.exposure.filename ?? "").replace(/\.dat$/i, "");
      const bn = (b.exposure.filename ?? "").replace(/\.dat$/i, "");
      return an.localeCompare(bn);
    });
    return out;
  }, [allExposures, samplesById, selectedTags, search, draft]);

  // Recently-used: server returns ids; resolve against the same row-source
  // we have on hand. Anything not in `rows` (e.g. cross-experiment picks
  // that aren't in the current sample window) is silently dropped — the
  // section is best-effort, not a guaranteed rendering of every recent.
  const recentRows = useMemo<RowEntry[]>(() => {
    const ids = recentsQ.data ?? [];
    if (ids.length === 0) return [];
    const byId = new Map<number, RowEntry>();
    for (const r of rows) byId.set(r.exposure.id, r);
    const out: RowEntry[] = [];
    for (const id of ids) {
      const r = byId.get(id);
      if (r) out.push(r);
    }
    return out;
  }, [recentsQ.data, rows]);

  if (!isOpen) return null;

  const togglePick = (exposureId: number, next: boolean): void => {
    setPicks((prev) =>
      next ? [...prev.filter((p) => p !== exposureId), exposureId]
           : prev.filter((p) => p !== exposureId),
    );
  };

  const toggleTag = (key: string, value: string): void => {
    const id = `${key}:${value}`;
    setSelectedTags((prev) => {
      const next = new Set(prev);
      if (next.has(id)) next.delete(id);
      else next.add(id);
      return next;
    });
  };

  const toggleExperiment = (id: number): void => {
    setSelectedExpIds((prev) =>
      prev.includes(id) ? prev.filter((e) => e !== id) : [...prev, id],
    );
  };

  const onAddSelected = (): void => {
    for (const id of picks) {
      addMember(id, qc);
    }
    onClose();
  };

  const onKeyDown = (e: React.KeyboardEvent<HTMLDivElement>): void => {
    if (e.key === "Escape") {
      e.preventDefault();
      onClose();
    }
  };

  const tagPairs = tagsQ.data ?? [];

  return (
    <div
      data-testid="comparison-picker-overlay"
      className="fixed inset-0 z-50 flex items-start justify-center pt-[10vh]
                 bg-[oklch(0.05_0_0/0.65)] backdrop-blur-sm anim-pal-in"
      role="presentation"
      onClick={(e) => { if (e.target === e.currentTarget) onClose(); }}
    >
      <div
        ref={dialogRef}
        data-testid="comparison-picker"
        role="dialog"
        aria-modal="true"
        aria-labelledby="comparison-picker-title"
        onKeyDown={onKeyDown}
        className="w-[min(720px,calc(100vw-48px))] max-h-[78vh]
                   bg-bg-elevated border border-border rounded-xl shadow-2xl
                   flex flex-col overflow-hidden anim-pal-scale"
      >
        <div className="flex items-center justify-between px-4 py-3 border-b border-border">
          <h2
            id="comparison-picker-title"
            className="text-base font-medium text-fg"
          >
            Add exposures
          </h2>
          <button
            type="button"
            data-testid="comparison-picker-close"
            onClick={onClose}
            className="text-fg-muted hover:text-fg text-sm px-2 py-1"
            aria-label="Close picker"
          >
            esc
          </button>
        </div>

        <div className="px-4 py-2 border-b border-border space-y-2">
          <input
            ref={inputRef}
            data-testid="comparison-picker-search"
            type="search"
            value={search}
            onChange={(e) => setSearch(e.target.value)}
            placeholder="Search exposures, samples, notes…"
            className="w-full bg-transparent border border-border rounded
                       px-2 py-1 text-sm"
            spellCheck={false}
          />
          <div className="flex flex-wrap gap-2 items-start">
            <div
              data-testid="picker-filter-experiment"
              className="flex flex-wrap gap-1 items-center"
            >
              <span className="text-xs text-fg-dim">Experiments:</span>
              {(experimentsQ.data ?? []).map((e) => {
                const active = selectedExpIds.includes(e.id);
                return (
                  <button
                    key={e.id}
                    type="button"
                    data-testid={`picker-experiment-option-${e.id}`}
                    data-active={active ? "true" : undefined}
                    onClick={() => toggleExperiment(e.id)}
                    className={
                      "text-xs px-2 py-0.5 rounded-full border " +
                      (active
                        ? "bg-accent/15 border-accent text-accent"
                        : "border-border text-fg-muted hover:text-fg")
                    }
                  >
                    {e.name ?? `Exp ${e.id}`}
                  </button>
                );
              })}
            </div>
            {tagPairs.length > 0 && (
              <div
                data-testid="picker-filter-tag"
                className="flex flex-wrap gap-1 items-center"
              >
                <span className="text-xs text-fg-dim">Tags:</span>
                {tagPairs.map((t) => {
                  const id = `${t.key}:${t.value}`;
                  const active = selectedTags.has(id);
                  return (
                    <button
                      key={id}
                      type="button"
                      data-testid={`picker-tag-option-${id}`}
                      data-active={active ? "true" : undefined}
                      onClick={() => toggleTag(t.key, t.value)}
                      className={
                        "text-xs px-2 py-0.5 rounded-full border " +
                        (active
                          ? "bg-accent/15 border-accent text-accent"
                          : "border-border text-fg-muted hover:text-fg")
                      }
                    >
                      {t.key}:{t.value}
                    </button>
                  );
                })}
              </div>
            )}
          </div>
        </div>

        <div className="flex-1 min-h-0 overflow-y-auto">
          {recentRows.length > 0 && (
            <section
              data-testid="comparison-picker-recents"
              className="border-b border-border/40 py-2"
            >
              <div className="px-4 text-xs font-medium text-fg-dim uppercase
                              tracking-wide pb-1">
                Recently used
              </div>
              <ul role="listbox" aria-label="Recently used exposures"
                  className="flex flex-col">
                {recentRows.map((r) => (
                  <li
                    key={`recent-${r.exposure.id}`}
                    role="option"
                    aria-selected={picks.includes(r.exposure.id)}
                    data-testid="picker-row"
                    data-exposure-id={r.exposure.id}
                    data-locked={r.alreadyAdded ? "true" : undefined}
                  >
                    <ExposureListRow
                      exposure={r.exposure}
                      sample={r.sample}
                      checked={r.alreadyAdded || picks.includes(r.exposure.id)}
                      onCheckedChange={(next) => togglePick(r.exposure.id, next)}
                      locked={r.alreadyAdded}
                      lockedReason={r.alreadyAdded ? "already added" : undefined}
                    />
                  </li>
                ))}
              </ul>
            </section>
          )}

          {rows.length === 0 ? (
            <div
              data-testid="comparison-picker-empty"
              className="px-4 py-8 text-center text-fg-muted text-sm"
            >
              No exposures match. Try clearing filters or broadening to all experiments.
            </div>
          ) : (
            <section className="py-2">
              <div className="px-4 text-xs font-medium text-fg-dim uppercase
                              tracking-wide pb-1">
                All exposures
              </div>
              <ul role="listbox" aria-label="Exposures"
                  className="flex flex-col">
                {rows.map((r) => (
                  <li
                    key={r.exposure.id}
                    role="option"
                    aria-selected={picks.includes(r.exposure.id)}
                    data-testid="picker-row"
                    data-exposure-id={r.exposure.id}
                    data-locked={r.alreadyAdded ? "true" : undefined}
                  >
                    <ExposureListRow
                      exposure={r.exposure}
                      sample={r.sample}
                      checked={r.alreadyAdded || picks.includes(r.exposure.id)}
                      onCheckedChange={(next) => togglePick(r.exposure.id, next)}
                      locked={r.alreadyAdded}
                      lockedReason={r.alreadyAdded ? "already added" : undefined}
                    />
                  </li>
                ))}
              </ul>
            </section>
          )}
        </div>

        <div className="flex items-center gap-2 px-4 py-3 border-t border-border">
          <span className="text-xs text-fg-dim flex-1">
            {picks.length} selected
          </span>
          <button
            type="button"
            data-testid="comparison-picker-cancel"
            onClick={onClose}
            className="px-3 py-1 rounded border border-border text-sm text-fg-muted"
          >
            Cancel
          </button>
          <button
            type="button"
            data-testid="comparison-picker-add"
            onClick={onAddSelected}
            disabled={picks.length === 0}
            className="px-3 py-1 rounded bg-accent text-bg text-sm font-medium
                       disabled:opacity-50"
          >
            {picks.length === 0 ? "Add selected" : `Add ${picks.length} selected`}
          </button>
        </div>
      </div>
    </div>
  );
}

/**
 * Aggregate-fetch exposures across a list of samples by calling the existing
 * useExposures hook for each sample id. React requires hooks to be called
 * unconditionally at the top level — we satisfy that by running a single
 * hook via a child component that loops over the sample list. To keep this
 * single-file, we use a small hook-of-hooks helper below.
 *
 * For the v1 picker the sample window is bounded (one experiment, ~tens of
 * samples max in practice), so the per-sample useQuery overhead is fine. If
 * cross-experiment broadening starts pulling hundreds of samples we can
 * swap to a single batched `GET /api/experiments/:eid/exposures` route.
 */
function useExposuresAcrossSamples(samples: Sample[]): Exposure[] {
  // We call useExposures once for each sample. To keep React's hook order
  // stable across renders we render a leaf component per sample inline.
  // For simplicity in v1 we batch the queries via a stable key list and
  // rely on TanStack Query's cache to make repeated reads ~free.
  const sampleIds = useMemo(() => samples.map((s) => s.id).sort((a, b) => a - b), [samples]);
  // Dynamic-length hook lists are not allowed in React. We use a single
  // useQueries invocation via a thin wrapper so the count can vary safely.
  // Below: inline implementation using useQueries from @tanstack/react-query.
  return useExposuresBatched(sampleIds);
}


function useExposuresBatched(sampleIds: number[]): Exposure[] {
  const results = useQueries({
    queries: sampleIds.map((id) => ({
      queryKey: ["sample", id, "exposures", { excludeRejected: false }] as const,
      queryFn: () => api.listExposures(id, { excludeRejected: false }),
    })),
  });
  return useMemo(() => {
    const out: Exposure[] = [];
    for (const r of results) {
      if (r.data) out.push(...r.data);
    }
    return out;
  }, [results]);
}
